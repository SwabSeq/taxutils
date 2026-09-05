use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::mpsc::{self, RecvTimeoutError};
use std::sync::Mutex;
use std::thread;
use std::time::Duration;

use anyhow::{Error, Result as AnyResult};
use pyo3::create_exception;
use pyo3::exceptions::{PyException, PyOSError, PyValueError};
use pyo3::prelude::*;
use taxutils_core::{self as core, CancellationToken, FilterMode};

const BACKEND_API_VERSION: u32 = 2;

create_exception!(_rust, TaxutilsBackendError, PyException);

fn python_error(error: Error) -> PyErr {
    if let Some(io_error) = error.downcast_ref::<std::io::Error>() {
        return PyOSError::new_err(io_error.to_string());
    }

    let message = format!("{error:#}");
    if message.starts_with("--batch-size")
        || message.starts_with("No accession found")
        || message.starts_with("No accessions were found")
        || message.starts_with("Invalid taxid")
        || message.contains("did not contain any taxids")
    {
        PyValueError::new_err(message)
    } else {
        TaxutilsBackendError::new_err(message)
    }
}

fn interruptible<T, F>(py: Python<'_>, operation: F) -> PyResult<T>
where
    T: Send + 'static,
    F: FnOnce(CancellationToken) -> AnyResult<T> + Send + 'static,
{
    let cancellation = CancellationToken::default();
    let worker_cancellation = cancellation.clone();
    let (sender, receiver) = mpsc::sync_channel(1);
    let worker = thread::spawn(move || {
        let _ = sender.send(operation(worker_cancellation));
    });
    let receiver = Mutex::new(receiver);
    let mut interrupt = None;

    let result = loop {
        match py.allow_threads(|| {
            receiver
                .lock()
                .expect("backend result receiver lock is not poisoned")
                .recv_timeout(Duration::from_millis(50))
        }) {
            Ok(result) => break result,
            Err(RecvTimeoutError::Timeout) => {
                if interrupt.is_none() {
                    if let Err(error) = py.check_signals() {
                        cancellation.cancel();
                        interrupt = Some(error);
                    }
                }
            }
            Err(RecvTimeoutError::Disconnected) => {
                let _ = worker.join();
                return Err(TaxutilsBackendError::new_err(
                    "the Rust backend worker stopped unexpectedly",
                ));
            }
        }
    };
    let _ = worker.join();
    if let Some(error) = interrupt {
        return Err(error);
    }
    result.map_err(python_error)
}

#[pyfunction]
fn api_version() -> u32 {
    BACKEND_API_VERSION
}

#[pyfunction]
fn crate_version() -> &'static str {
    core::VERSION
}

#[pyfunction]
#[pyo3(signature = (fasta_path, output_path, batch_size=10_000))]
fn extract_accessions(
    py: Python<'_>,
    fasta_path: PathBuf,
    output_path: PathBuf,
    batch_size: usize,
) -> PyResult<usize> {
    interruptible(py, move |cancellation| {
        core::extract_accessions_with_cancel(fasta_path, output_path, batch_size, &cancellation)
    })
}

#[pyfunction]
#[pyo3(signature = (input_path, output_path=None, verbose=false))]
fn clean_fasta_headers(
    py: Python<'_>,
    input_path: PathBuf,
    output_path: Option<PathBuf>,
    verbose: bool,
) -> PyResult<()> {
    interruptible(py, move |cancellation| {
        core::clean_fasta_headers_with_cancel(
            input_path,
            output_path.as_deref(),
            verbose,
            &cancellation,
        )
    })
}

#[pyfunction]
#[pyo3(signature = (input_path, accession_query, output_path, version=true, batch_size=1_000_000, verbose=false))]
fn grep_fasta(
    py: Python<'_>,
    input_path: PathBuf,
    accession_query: String,
    output_path: PathBuf,
    version: bool,
    batch_size: usize,
    verbose: bool,
) -> PyResult<(usize, usize, usize, usize)> {
    let stats = interruptible(py, move |cancellation| {
        core::grep_fasta_with_cancel(
            input_path,
            &accession_query,
            output_path,
            version,
            batch_size,
            verbose,
            &cancellation,
        )
    })?;
    Ok((
        stats.requested,
        stats.scanned,
        stats.matched,
        stats.missing_accession,
    ))
}

#[pyfunction]
#[pyo3(signature = (input_path, output_path, filter_taxa, save_folder, filter_mode="remove", batch_size=5_000, verbose=false, wgs=false))]
#[allow(clippy::too_many_arguments)]
fn filter_fasta(
    py: Python<'_>,
    input_path: PathBuf,
    output_path: Option<PathBuf>,
    filter_taxa: Vec<i64>,
    save_folder: PathBuf,
    filter_mode: &str,
    batch_size: usize,
    verbose: bool,
    wgs: bool,
) -> PyResult<(usize, usize, usize, usize)> {
    let mode = match filter_mode {
        "keep" => FilterMode::Keep,
        "remove" => FilterMode::Remove,
        _ => {
            return Err(PyValueError::new_err(
                "filter_mode must be 'keep' or 'remove'",
            ))
        }
    };
    let taxa = filter_taxa.into_iter().collect::<HashSet<_>>();
    let stats = interruptible(py, move |cancellation| {
        core::filter_fasta_with_options_and_cancel(
            input_path,
            output_path.as_deref(),
            &taxa,
            mode,
            batch_size,
            verbose,
            save_folder,
            wgs,
            &cancellation,
        )
    })?;
    Ok((
        stats.kept,
        stats.removed,
        stats.missing_accession,
        stats.missing_taxid,
    ))
}

#[pyfunction]
#[pyo3(signature = (save_folder, rebuild=false, wgs=false, keep_downloads=true))]
fn ensure_accession_database(
    py: Python<'_>,
    save_folder: PathBuf,
    rebuild: bool,
    wgs: bool,
    keep_downloads: bool,
) -> PyResult<PathBuf> {
    py.allow_threads(move || {
        core::ensure_accession_database(
            save_folder,
            core::AccessionDatabaseOptions {
                rebuild,
                wgs,
                keep_downloads,
            },
        )
    })
    .map_err(python_error)
}

#[pyfunction]
#[pyo3(signature = (save_folder, accessions, low_memory, wgs, keep_downloads=true))]
fn lookup_accession_taxids(
    py: Python<'_>,
    save_folder: PathBuf,
    accessions: Vec<String>,
    low_memory: bool,
    wgs: bool,
    keep_downloads: bool,
) -> PyResult<HashMap<String, i64>> {
    let requested = accessions.into_iter().collect::<HashSet<_>>();
    if requested.is_empty() {
        return Ok(HashMap::new());
    }
    py.allow_threads(move || {
        if !low_memory {
            core::ensure_accession_database(
                &save_folder,
                core::AccessionDatabaseOptions {
                    wgs,
                    keep_downloads,
                    ..Default::default()
                },
            )?;
        }
        core::lookup_accession_taxids(save_folder, requested, low_memory, wgs)
    })
    .map_err(python_error)
}

#[pyfunction]
#[pyo3(signature = (save_folder, taxa, low_memory, wgs, keep_downloads=true))]
fn lookup_taxid_accessions(
    py: Python<'_>,
    save_folder: PathBuf,
    taxa: Vec<i64>,
    low_memory: bool,
    wgs: bool,
    keep_downloads: bool,
) -> PyResult<HashSet<String>> {
    if taxa.is_empty() {
        return Ok(HashSet::new());
    }
    py.allow_threads(move || {
        if !low_memory {
            core::ensure_accession_database(
                &save_folder,
                core::AccessionDatabaseOptions {
                    wgs,
                    keep_downloads,
                    ..Default::default()
                },
            )?;
        }
        core::lookup_taxid_accessions(save_folder, &taxa, low_memory, wgs)
    })
    .map_err(python_error)
}

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "TaxutilsBackendError",
        module.py().get_type::<TaxutilsBackendError>(),
    )?;
    module.add_function(wrap_pyfunction!(api_version, module)?)?;
    module.add_function(wrap_pyfunction!(crate_version, module)?)?;
    module.add_function(wrap_pyfunction!(extract_accessions, module)?)?;
    module.add_function(wrap_pyfunction!(clean_fasta_headers, module)?)?;
    module.add_function(wrap_pyfunction!(grep_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(filter_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(ensure_accession_database, module)?)?;
    module.add_function(wrap_pyfunction!(lookup_accession_taxids, module)?)?;
    module.add_function(wrap_pyfunction!(lookup_taxid_accessions, module)?)?;
    Ok(())
}

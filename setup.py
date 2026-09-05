"""Build configuration for the required Rust backend."""

from setuptools import setup
from setuptools_rust import Binding, RustExtension


setup(
    rust_extensions=[
        RustExtension(
            "taxutils._rust",
            path="rust/Cargo.toml",
            binding=Binding.PyO3,
            py_limited_api=True,
        )
    ],
    zip_safe=False,
    options={"bdist_wheel": {"py_limited_api": "cp310"}},
)

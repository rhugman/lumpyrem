"""Reference-oracle fixtures for the lumpyrem port.

This package drives the original LUMPREM Fortran so that every phase of the
port can be gated on agreement with it.  It is the only place legacy file
formats are understood, and it is never imported by ``lumpyrem`` itself.
"""

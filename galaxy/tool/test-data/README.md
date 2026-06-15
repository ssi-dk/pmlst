# Galaxy Test Data

This directory contains an offline pMLST fixture for Galaxy wrapper tests.

- `test.fsa` is the existing repository FASTA fixture copied for Galaxy tests.
- `pmlst_fixture_db/` is an IncF-only subset of the full pMLST database.
- The fixture database includes root metadata plus `incf.fsa`,
  `incf.txt.clean`, `incf.clpx`, and `incf.name`.
- KMA indexes are intentionally omitted because Galaxy tests use FASTA/BLAST
  only.

The fixture is under 1 MB and must not be used as a production pMLST database.

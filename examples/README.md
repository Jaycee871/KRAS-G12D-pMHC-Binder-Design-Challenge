# examples/

These files exist **only to demonstrate the file format and to exercise the
tooling in CI**. They are *not* binder designs, they have not been modelled
against any target, their score columns are all zero, and they must never be
submitted to the competition. Every sequence here is named
`EXAMPLE_FORMAT_ONLY_*` for that reason.

* `scored_designs_example.csv` - the shape of the table `make_submission.py`
  consumes (`name`, `sequence`, `molecule_class`, plus score columns).
  The score columns are the gate inputs listed in `docs/PIPELINE.md`
  section 6 and checked by `src/design_gates.py`, including the decoy
  counter-screen (`iptm_decoy`) and the buried-Asp term (`asp6_buried_frac`).
* `example_submission.csv` - the shape of the file uploaded to Proteinbase.

Try them:

```bash
python src/make_submission.py examples/scored_designs_example.csv -o /tmp/sub.csv -n 3
python src/validate_submission.py /tmp/sub.csv
```

## Environmnet Setup

```bash
conda create -y -n yarp python=3.12
conda activate yarp
pip install rdkit scipy
```

## Test Run

```bash
python enumerate.py -s "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO" -m b2f2 -o test.json
```

test.json should look like:
```json
[
  "O=CC(O)C(O)C(O)C(O)CO>>O.O=CC(O)=CC(O)C(O)CO",
  ...
]
```

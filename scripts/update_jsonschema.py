"""Update `schema.json` with latest output schema.

% python3 ./scripts/update_jsonschema.py
"""

import json
from pathlib import Path

from dcd_mapping.schemas import ScoresetMapping

schema = ScoresetMapping.model_json_schema()
schema_text = json.dumps(schema, indent=4)
file = Path(__file__).parents[1] / "schema.json"

with file.open("w") as f:
    f.write(schema_text)

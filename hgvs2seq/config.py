from pydantic import BaseModel
import yaml
from typing import Optional, List, Tuple

class TranscriptConfig(BaseModel):
    """
    Configuration for a single transcript, including its structure and genomic context.
    """
    transcript_id: str
    gene_symbol: Optional[str] = None
    assembly: str
    strand: int
    exons: List[Tuple[int, int]]
    cds_start_c: Optional[int] = None
    cds_end_c: Optional[int] = None

def load_config(path: str) -> TranscriptConfig:
    """
    Loads transcript configuration from a JSON or YAML file.
    """
    with open(path, 'r') as f:
        if path.endswith(".json"):
            import json
            data = json.load(f)
        elif path.endswith((".yaml", ".yml")):
            data = yaml.safe_load(f)
        else:
            raise ValueError("Configuration file must be a .json or .yaml file")
    return TranscriptConfig(**data)
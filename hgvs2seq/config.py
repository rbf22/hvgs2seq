"""
Handles loading and validation of the transcript configuration.
"""
import json
import logging
from .models import TranscriptConfig
from pydantic import ValidationError

_logger = logging.getLogger(__name__)

def load_config(config_path: str) -> TranscriptConfig:
    """
    Loads a transcript configuration from a JSON file and validates it.

    Args:
        config_path: The path to the JSON configuration file.

    Returns:
        A validated TranscriptConfig object.

    Raises:
        FileNotFoundError: If the config file does not exist.
        ValueError: If the config file is not valid JSON or fails validation.
    """
    _logger.info(f"Loading transcript configuration from {config_path}")
    try:
        with open(config_path, 'r') as f:
            config_data = json.load(f)
    except FileNotFoundError:
        _logger.error(f"Configuration file not found at {config_path}")
        raise
    except json.JSONDecodeError as e:
        _logger.error(f"Invalid JSON in configuration file: {config_path}")
        raise ValueError(f"Invalid JSON in {config_path}: {e}") from e

    try:
        config = TranscriptConfig(**config_data)
        _logger.info(f"Successfully loaded and validated configuration for {config.transcript_id}")
        return config
    except ValidationError as e:
        _logger.error("Configuration validation failed.")
        raise ValueError(f"Configuration validation failed: {e}") from e
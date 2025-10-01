"""Pytest configuration and fixtures for tests."""
import pytest
from dataclasses import dataclass

@dataclass
class MockTranscriptConfig:
    """Mock transcript configuration for testing."""
    transcript_id: str = "test_transcript"
    cds_start: int = 1
    cds_end: int = 1000
    exons: list = None
    
    def __post_init__(self):
        if self.exons is None:
            self.exons = [(1, 1000)]

@pytest.fixture
def create_mock_config():
    """Create a mock transcript configuration."""
    def _create_mock_config(transcript_id="test_transcript", cds_start=1, cds_end=1000, exons=None):
        return MockTranscriptConfig(
            transcript_id=transcript_id,
            cds_start=cds_start,
            cds_end=cds_end,
            exons=exons or [(1, 1000)]
        )
    return _create_mock_config

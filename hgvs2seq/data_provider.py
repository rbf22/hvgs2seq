"""
Initializes and holds connections to data sources (UTA, SeqRepo).
This avoids reconnecting repeatedly and makes providers accessible across the package.
Uses lazy initialization to avoid connecting at import time, which is crucial for testing.
"""
import hgvs.dataproviders.uta
import hgvs.assemblymapper
from biocommons.seqrepo import SeqRepo
import os
import logging

_logger = logging.getLogger(__name__)

# Global variables to hold the initialized connections (singletons)
_hdp = None
_sr = None
_am = None

def get_hdp():
    """Returns a singleton HGVS data provider (UTA) instance."""
    global _hdp
    if _hdp is None:
        _logger.info("Initializing connection to UTA database...")
        try:
            _hdp = hgvs.dataproviders.uta.connect()
            _logger.info("Successfully connected to UTA database.")
        except Exception as e:
            _logger.error("Failed to connect to UTA database.")
            raise ConnectionError(
                "Failed to connect to UTA database. Ensure that the HGVS_UTA_URL "
                "environment variable is set correctly and the database is accessible."
            ) from e
    return _hdp

def get_sr():
    """Returns a singleton SeqRepo instance."""
    global _sr
    if _sr is None:
        _logger.info("Initializing connection to SeqRepo...")
        seqrepo_root = os.environ.get("SEQREPO_ROOT")
        if not seqrepo_root or not os.path.exists(seqrepo_root):
            _logger.error("SEQREPO_ROOT environment variable not set or path does not exist.")
            raise EnvironmentError(
                "SEQREPO_ROOT environment variable is not set or path is invalid. "
                "Please set it to the path of your seqrepo instance."
            )
        _sr = SeqRepo(seqrepo_root)
        _logger.info(f"Successfully connected to SeqRepo at {seqrepo_root}.")
    return _sr

def get_am():
    """Returns a singleton AssemblyMapper instance."""
    global _am
    if _am is None:
        _logger.info("Initializing AssemblyMapper...")
        _am = hgvs.assemblymapper.AssemblyMapper(
            get_hdp(),
            seqrepo=get_sr(),
            replace_reference=True
        )
        _logger.info("AssemblyMapper initialized.")
    return _am
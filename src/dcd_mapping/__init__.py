"""Provide tools for mapping MaveDB experiments to human reference sequences for
use in scientific and clinical applications. A project of the Data Coordination and
Dissemination (DCD) workstream of the Atlas of Variant Effects (AVE) Alliance.

See the mapping manuscript for more information:
https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1
"""

from dotenv import load_dotenv

from .main import map_scoreset, map_scoreset_urn

__all__ = ["map_scoreset", "map_scoreset_urn"]

load_dotenv()

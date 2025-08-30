"""
Servicio directo de PubMed usando Entrez API
"""

import sys
import os
from mcp.server.fastmcp import FastMCP
from typing import List, Dict, Any
from .pubmed_service import PubMedService

# Agregar el directorio padre al path para importar configuraciones
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))



mcp = FastMCP(__name__, 'pubmed')
pubmed_service = PubMedService()


@mcp.tool()
async def search_pubmed(query: str) -> List[Dict[str, Any]]:
    """Busca artículos en PubMed"""
    return await pubmed_service.search_and_fetch(query)
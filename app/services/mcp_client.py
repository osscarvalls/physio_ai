"""
Cliente MCP para conectar con el servidor de PubMed
"""

import asyncio
import logging
from typing import List, Dict, Any, Optional
from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters
from mcp.types import TextContent

logger = logging.getLogger(__name__)


class PubMedMCPClient:
    """Cliente MCP para conectar con el servidor de PubMed"""
    
    def __init__(self, server_path: str = None):
        """
        Inicializa el cliente MCP
        
        Args:
            server_path: Ruta al ejecutable del servidor MCP (opcional)
        """
        self.server_path = server_path or "python -m pubmed_mcp.server"
        self.session: Optional[ClientSessionGroup] = None
    
    async def connect(self):
        """Conecta con el servidor MCP"""
        try:
            # Configurar parámetros del servidor
            server_params = StdioServerParameters(
                command=self.server_path.split(),
                env={}
            )
            
            # Crear sesión del cliente
            self.session = ClientSessionGroup(server_params)
            
            # Inicializar la sesión
            await self.session.initialize()
            
            logger.info("✅ Cliente MCP conectado exitosamente")
            
        except Exception as e:
            logger.error(f"❌ Error conectando con el servidor MCP: {e}")
            raise
    
    async def disconnect(self):
        """Desconecta del servidor MCP"""
        if self.session:
            try:
                await self.session.close()
                logger.info("✅ Cliente MCP desconectado")
            except Exception as e:
                logger.error(f"❌ Error desconectando del servidor MCP: {e}")
            finally:
                self.session = None
    
    async def search_pubmed(self, query: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Busca artículos en PubMed usando el servidor MCP
        
        Args:
            query: Consulta de búsqueda
            max_results: Número máximo de resultados
            
        Returns:
            Lista de artículos encontrados
        """
        if not self.session:
            raise RuntimeError("Cliente MCP no conectado. Llama a connect() primero.")
        
        try:
            # Obtener herramientas disponibles
            tools = await self.session.list_tools()
            
            # Buscar la herramienta search_pubmed
            search_tool = None
            for tool in tools:
                if tool.name == 'search_pubmed':
                    search_tool = tool
                    break
            
            if not search_tool:
                raise RuntimeError("Herramienta 'search_pubmed' no encontrada en el servidor MCP")
            
            # Ejecutar la búsqueda
            result = await self.session.call_tool(search_tool.name, {"query": query})
            
            # Procesar la respuesta
            if result and hasattr(result, 'articles'):
                articles = result.articles
            elif isinstance(result, dict) and 'articles' in result:
                articles = result['articles']
            elif isinstance(result, (list, tuple)) and len(result) > 0:
                articles = result[0] if isinstance(result[0], list) else result
            else:
                articles = []
            
            # Limitar resultados si es necesario
            if max_results and len(articles) > max_results:
                articles = articles[:max_results]
            
            logger.info(f"✅ Búsqueda MCP completada: {len(articles)} artículos encontrados")
            return articles
            
        except Exception as e:
            logger.error(f"❌ Error en búsqueda MCP: {e}")
            raise
    
    async def __aenter__(self):
        """Context manager para conectar automáticamente"""
        await self.connect()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Context manager para desconectar automáticamente"""
        await self.disconnect()

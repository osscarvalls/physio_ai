"""
Servicio de búsqueda semántica para evidencia médica usando langchain-qdrant
"""

import logging
from typing import List, Dict, Any
from langchain_openai import OpenAIEmbeddings
from langchain_qdrant import QdrantVectorStore
from qdrant_client import QdrantClient
from qdrant_client.models import Distance, VectorParams

from app.config.settings import settings

logger = logging.getLogger(__name__)


class SemanticSearchService:
    """Servicio de búsqueda semántica para evidencia médica usando langchain-qdrant"""
    
    def __init__(self):
        """Inicializa el servicio de búsqueda semántica"""
        self._setup_components()
    
    def _setup_components(self):
        """Configura los componentes de LangChain con Qdrant"""
        try:
            # Configurar embeddings
            self.embeddings = OpenAIEmbeddings(
                api_key=settings.OPENAI_API_KEY
            )
            
            # Crear cliente de Qdrant para verificar/crear colección
            self.qdrant_client = QdrantClient(
                url=settings.QDRANT_URL,
                api_key=settings.QDRANT_API_KEY
            )
            
            # Verificar si la colección existe, si no, crearla
            self._ensure_collection_exists()
            
            # Configurar vector store de Qdrant
            self.vector_store = QdrantVectorStore.from_existing_collection(
                collection_name=settings.QDRANT_COLLECTION_NAME,
                embedding=self.embeddings,
                url=settings.QDRANT_URL,
                api_key=settings.QDRANT_API_KEY
            )
            
            logger.info("Componentes de búsqueda semántica configurados correctamente")
            
        except Exception as e:
            logger.error(f"Error configurando componentes: {str(e)}")
            self.embeddings = None
            self.vector_store = None
            self.qdrant_client = None
    
    def _ensure_collection_exists(self):
        """Asegura que la colección existe, creándola si es necesario"""
        try:
            # Verificar si la colección existe
            collections = self.qdrant_client.get_collections()
            collection_names = [col.name for col in collections.collections]
            
            if settings.QDRANT_COLLECTION_NAME not in collection_names:
                logger.info(f"Creando colección '{settings.QDRANT_COLLECTION_NAME}' en Qdrant")
                
                # Crear la colección con configuración estándar para embeddings
                self.qdrant_client.create_collection(
                    collection_name=settings.QDRANT_COLLECTION_NAME,
                    vectors_config=VectorParams(
                        size=1536,  # Tamaño estándar para OpenAI embeddings
                        distance=Distance.COSINE
                    )
                )
                
                logger.info(f"Colección '{settings.QDRANT_COLLECTION_NAME}' creada exitosamente")
            else:
                logger.info(f"Colección '{settings.QDRANT_COLLECTION_NAME}' ya existe")
                
        except Exception as e:
            logger.error(f"Error asegurando que la colección existe: {str(e)}")
            raise
    
    async def search_medical_evidence(self, symptoms: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Busca evidencia médica relevante para los síntomas usando Qdrant
        
        Args:
            symptoms: Síntomas del paciente
            max_results: Número máximo de resultados
            
        Returns:
            Lista de documentos relevantes
        """
        try:
            if not self.vector_store:
                logger.warning("Vector store no disponible")
                return []
            
            # Verificar si la colección tiene documentos
            collection_info = self.qdrant_client.get_collection(settings.QDRANT_COLLECTION_NAME)
            if collection_info.points_count == 0:
                logger.info("La colección está vacía, no hay evidencia médica disponible")
                return []
            
            # Búsqueda semántica en Qdrant
            results = self.vector_store.similarity_search(
                symptoms, 
                k=max_results
            )
            
            # Formatear resultados
            formatted_results = []
            for i, doc in enumerate(results):
                formatted_results.append({
                    'title': doc.metadata.get('title', 'Sin título'),
                    'content': doc.page_content,
                    'metadata': doc.metadata,
                    'rank': i + 1
                })
            
            logger.info(f"Búsqueda semántica completada, {len(formatted_results)} resultados")
            return formatted_results
            
        except Exception as e:
            logger.error(f"Error en búsqueda semántica: {str(e)}")
            return []

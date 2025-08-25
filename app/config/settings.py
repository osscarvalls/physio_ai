"""
Configuración centralizada de la aplicación PySIO AI
"""

import os
from typing import Optional
from pydantic import Field
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Configuración de la aplicación PySIO AI"""
    
    # Configuración del servidor
    HOST: str = Field(default="0.0.0.0", description="Host del servidor")
    PORT: int = Field(default=8000, description="Puerto del servidor")
    DEBUG: bool = Field(default=False, description="Modo debug")
    
    # Configuración de OpenAI
    OPENAI_API_KEY: str = Field(..., description="API Key de OpenAI")
    OPENAI_PROJECT_ID: str = Field(..., description="ID de proyecto en OPENAI")
    
    # Configuración de PubMed
    ENTREZ_EMAIL: str = Field(default="user@example.com", description="Email para Entrez/PubMed")
    
    # Configuración de QDRANT CLOUD
    QDRANT_URL: str = Field(..., description="URL de Qdrant Cloud (ej: https://your-cluster.qdrant.io)")
    QDRANT_API_KEY: str = Field(..., description="API Key de Qdrant Cloud")
    QDRANT_COLLECTION_NAME: str = Field(default="pysio_ai", description="Nombre de la colección en Qdrant")
    
    # Configuración de logging
    LOG_LEVEL: str = Field(default="DEBUG", description="Nivel de logging")
    
    class Config:
        env_file = ".env"
        case_sensitive = True


# Instancia global de configuración
settings = Settings()

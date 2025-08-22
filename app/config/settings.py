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
    OPENAI_ORGANIZATION_ID: Optional[str] = Field(default=None, description="ID de organización de OpenAI")
    OPENAI_MODEL: str = Field(default="gpt-4o-mini", description="Modelo de OpenAI a usar")
    OPENAI_TEMPERATURE: float = Field(default=0.7, description="Temperatura para la generación")
    
    # Configuración de PubMed
    ENTREZ_EMAIL: str = Field(default="user@example.com", description="Email para Entrez/PubMed")
    
    # Configuración de Chroma DB
    CHROMA_PERSIST_DIR: str = Field(default="./chroma_db", description="Directorio de persistencia de Chroma DB")
    
    # Configuración de logging
    LOG_LEVEL: str = Field(default="INFO", description="Nivel de logging")
    
    class Config:
        env_file = ".env"
        case_sensitive = True


# Instancia global de configuración
settings = Settings()

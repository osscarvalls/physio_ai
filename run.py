#!/usr/bin/env python3
"""
Script simple para ejecutar la API PySIO AI
"""

import uvicorn
from app.config.settings import settings
from pubmed_mcp import mcp

if __name__ == "__main__":
    print("🚀 Iniciando PySIO AI API...")
    print(f"📊 Configuración: HOST={settings.HOST}, PORT={settings.PORT}")
    print(f"🔗 API Docs: http://{settings.HOST}:{settings.PORT}/docs")
    print(f"🔗 Health Check: http://{settings.HOST}:{settings.PORT}/api/v1/health")
    print("=" * 50)
    
        # MCP comentado temporalmente para evitar conflictos
    # print("ℹ️  MCP deshabilitado temporalmente")
    
    # Ejecutar la API FastAPI
    uvicorn.run(
        "app.main:app",
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.DEBUG,
        log_level=settings.LOG_LEVEL.lower()
    )
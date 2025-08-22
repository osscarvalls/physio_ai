#!/usr/bin/env python3
"""
Script de inicio para PySIO AI
"""

import sys
import os

# Agregar el directorio raíz al path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.main import app

if __name__ == "__main__":
    import uvicorn
    from app.config.settings import settings
    
    print("🏥 PySIO AI - Servicio de Diagnóstico Médico")
    print("=" * 50)
    print(f"📊 Configuración: {settings.HOST}:{settings.PORT}")
    print(f"🔧 Modo debug: {settings.DEBUG}")
    print(f"📚 Documentación: http://{settings.HOST}:{settings.PORT}/docs")
    print("=" * 50)
    
    uvicorn.run(
        "app.main:app",
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.DEBUG,
        log_level="info"
    )

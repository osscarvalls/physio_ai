"""
PySIO AI - Aplicación principal simplificada
"""

import logging
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.controllers.diagnosis_controller import diagnosis_router

# Configuración simple de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)

logger = logging.getLogger(__name__)

def create_app() -> FastAPI:
    """Crea y configura la aplicación FastAPI"""
    
    app = FastAPI(
        title="PySIO AI - API de Diagnóstico Médico",
        description="API inteligente para diagnóstico médico basado en síntomas usando IA",
        version="1.0.0"
    )
    
    # Configuración de CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    
    # Incluir routers
    app.include_router(diagnosis_router)
    
    # Endpoint raíz
    @app.get("/")
    async def root():
        return {
            "message": "PySIO AI - API de Diagnóstico Médico",
            "version": "1.0.0",
            "docs": "/docs"
        }
    
    return app


# Crear instancia de la aplicación
app = create_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)

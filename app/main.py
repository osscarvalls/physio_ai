"""
PySIO AI - Aplicación principal refactorizada
"""

import logging
import sys
from contextlib import asynccontextmanager

import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from app.config.settings import settings
from app.controllers.diagnosis_controller import diagnosis_router


# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("pysio_ai.log")
    ]
)

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Maneja el ciclo de vida de la aplicación"""
    # Startup
    logger.info("🚀 Iniciando PySIO AI...")
    logger.info(f"📊 Configuración: {settings.HOST}:{settings.PORT}")
    
    yield
    
    # Shutdown
    logger.info("🛑 Cerrando PySIO AI...")


def create_app() -> FastAPI:
    """Crea y configura la aplicación FastAPI"""
    
    app = FastAPI(
        title="PySIO AI",
        description="Asistente de diagnóstico médico basado en IA",
        version="1.0.0",
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=lifespan
    )
    
    # Configurar CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # En producción, especificar dominios específicos
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    
    # Montar archivos estáticos
    app.mount("/static", StaticFiles(directory=settings.STATIC_DIR), name="static")
    
    # Incluir routers
    app.include_router(diagnosis_router)
    
    # Endpoint raíz adicional
    @app.get("/")
    async def root():
        return {
            "message": "Bienvenido a PySIO AI",
            "version": "1.0.0",
            "docs": "/docs",
            "health": "/api/v1/health"
        }
    
    return app


# Crear instancia de la aplicación
app = create_app()


if __name__ == "__main__":
    logger.info("🏥 PySIO AI - Servicio de Diagnóstico Médico")
    logger.info("=" * 50)
    
    try:
        uvicorn.run(
            "app.main:app",
            host=settings.HOST,
            port=settings.PORT,
            reload=settings.DEBUG,
            log_level="info"
        )
    except KeyboardInterrupt:
        logger.info("Aplicación interrumpida por el usuario")
    except Exception as e:
        logger.error(f"Error iniciando la aplicación: {str(e)}")
        sys.exit(1)

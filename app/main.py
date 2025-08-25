"""
PySIO AI - Aplicación principal refactorizada
"""

import logging
import structlog
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.config.settings import settings
from app.controllers.diagnosis_controller import diagnosis_router

# Configuración de logging estructurado
structlog.configure(
    processors=[
        structlog.stdlib.filter_by_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
        structlog.processors.UnicodeDecoder(),
        structlog.processors.JSONRenderer()
    ],
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    wrapper_class=structlog.stdlib.BoundLogger,
    cache_logger_on_first_use=True,
)

logger = structlog.get_logger()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Maneja eventos de inicio y cierre de la aplicación"""
    # Startup
    logger.info("🚀 Iniciando PySIO AI API...")
    logger.info(f"📊 Configuración: HOST={settings.HOST}, PORT={settings.PORT}, DEBUG={settings.DEBUG}")
    
    yield
    
    # Shutdown
    logger.info("🛑 Cerrando PySIO AI API...")


def create_app() -> FastAPI:
    """Crea y configura la aplicación FastAPI"""
    
    app = FastAPI(
        title="PySIO AI - API de Diagnóstico Médico",
        description="API inteligente para diagnóstico médico basado en síntomas usando IA",
        version="1.0.0",
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=lifespan
    )
    
    # Configuración de CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # En producción, especificar dominios específicos
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
            "docs": "/docs",
            "health": "/api/v1/health"
        }
    
    return app


# Crear instancia de la aplicación
app = create_app()


if __name__ == "__main__":
    import uvicorn
    
    uvicorn.run(
        "app.main:app",
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.DEBUG,
        log_level=settings.LOG_LEVEL.lower()
    )

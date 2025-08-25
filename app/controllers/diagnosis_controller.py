"""
Controlador para endpoints de diagnóstico médico
"""

import logging
from fastapi import APIRouter, Depends, HTTPException
from typing import List, Dict, Any

from app.models.diagnosis import (
    DiagnosisRequest, 
    DiagnosisResponse, 
    MedicalEvidence
)
from app.services.diagnosis_service import DiagnosisService
from app.services.semantic_search_service import SemanticSearchService
from app.services.pubmed_service import PubMedService

logger = logging.getLogger(__name__)

# Router para los endpoints de diagnóstico
diagnosis_router = APIRouter(prefix="/api/v1", tags=["diagnosis"])

# Instancias singleton de los servicios
_diagnosis_service = None
_semantic_search_service = None
_pubmed_service = None


def get_diagnosis_service() -> DiagnosisService:
    """Dependency injection para DiagnosisService - Singleton pattern"""
    global _diagnosis_service
    if _diagnosis_service is None:
        logger.info("🔧 Configurando DiagnosisService por primera vez...")
        _diagnosis_service = DiagnosisService()
        logger.info("✅ DiagnosisService configurado y listo para uso")
    return _diagnosis_service


def get_semantic_search_service() -> SemanticSearchService:
    """Dependency injection para SemanticSearchService - Singleton pattern"""
    global _semantic_search_service
    if _semantic_search_service is None:
        logger.info("🔧 Configurando SemanticSearchService por primera vez...")
        _semantic_search_service = SemanticSearchService()
        logger.info("✅ SemanticSearchService configurado y listo para uso")
    return _semantic_search_service


def get_pubmed_service() -> PubMedService:
    """Dependency injection para PubMedService - Singleton pattern"""
    global _pubmed_service
    if _pubmed_service is None:
        logger.info("🔧 Configurando PubMedService por primera vez...")
        _pubmed_service = PubMedService()
        logger.info("✅ PubMedService configurado y listo para uso")
    return _pubmed_service


@diagnosis_router.post("/diagnosis", response_model=DiagnosisResponse)
async def create_diagnosis(
    request: DiagnosisRequest,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> DiagnosisResponse:
    """
    Crea un nuevo diagnóstico basado en síntomas y contexto médico
    """
    try:
        return await diagnosis_service.evaluate_patient(request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al generar diagnóstico: {str(e)}")

@diagnosis_router.get("/health")
async def health_check() -> Dict[str, str]:
    """
    Endpoint de verificación de salud de la API
    """
    return {"status": "healthy", "service": "PySIO AI Diagnosis API"}

@diagnosis_router.post("/reinitialize")
async def reinitialize_services() -> Dict[str, str]:
    """
    Reinicializa todos los servicios (útil para testing o cambios de configuración)
    """
    global _diagnosis_service, _semantic_search_service, _pubmed_service
    
    logger.info("🔄 Reinicializando todos los servicios...")
    
    # Limpiar instancias existentes
    _diagnosis_service = None
    _semantic_search_service = None
    _pubmed_service = None
    
    # Forzar nueva inicialización
    get_diagnosis_service()
    get_semantic_search_service()
    get_pubmed_service()
    
    logger.info("✅ Todos los servicios reinicializados")
    return {"status": "reinitialized", "message": "Todos los servicios han sido reinicializados"}

"""
Controlador para endpoints de diagnóstico médico
"""

import logging
from fastapi import APIRouter, Depends, HTTPException
from typing import List, Dict, Any

from app.models.diagnosis import DiagnosisRequest, DiagnosisResponse
from app.services.diagnosis_service import DiagnosisService
from app.services.semantic_search_service import SemanticSearchService

logger = logging.getLogger(__name__)

diagnosis_router = APIRouter(prefix="/api/v1", tags=["diagnosis"])

_diagnosis_service = None
_semantic_search_service = None
_pubmed_service = None


def get_diagnosis_service() -> DiagnosisService:
    """Dependency injection para DiagnosisService"""
    global _diagnosis_service
    if _diagnosis_service is None:
        logger.info("Configurando DiagnosisService...")
        _diagnosis_service = DiagnosisService()
        logger.info("DiagnosisService configurado")
    return _diagnosis_service


def get_semantic_search_service() -> SemanticSearchService:
    """Dependency injection para SemanticSearchService"""
    global _semantic_search_service
    if _semantic_search_service is None:
        logger.info("Configurando SemanticSearchService...")
        _semantic_search_service = SemanticSearchService()
        logger.info("SemanticSearchService configurado")
    return _semantic_search_service


def get_pubmed_service() -> PubMedService:
    """Dependency injection para PubMedService"""
    global _pubmed_service
    if _pubmed_service is None:
        logger.info("Configurando PubMedService...")
        _pubmed_service = PubMedService()
        logger.info("PubMedService configurado")
    return _pubmed_service


@diagnosis_router.post("/diagnosis", response_model=DiagnosisResponse)
async def create_diagnosis(
    request: DiagnosisRequest,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> DiagnosisResponse:
    """Crea un nuevo diagnóstico basado en síntomas y contexto médico"""
    try:
        return await diagnosis_service.evaluate_patient(request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al generar diagnóstico: {str(e)}")

@diagnosis_router.get("/health")
async def health_check() -> Dict[str, str]:
    """Endpoint de verificación de salud de la API"""
    return {"status": "healthy", "service": "PySIO AI Diagnosis API"}

@diagnosis_router.post("/reinitialize")
async def reinitialize_services() -> Dict[str, str]:
    """Reinicializa todos los servicios"""
    global _diagnosis_service, _semantic_search_service, _pubmed_service
    
    logger.info("Reinicializando servicios...")
    
    _diagnosis_service = None
    _semantic_search_service = None
    _pubmed_service = None
    
    get_diagnosis_service()
    get_semantic_search_service()
    get_pubmed_service()
    
    logger.info("Servicios reinicializados")
    return {"status": "reinitialized", "message": "Servicios reinicializados"}

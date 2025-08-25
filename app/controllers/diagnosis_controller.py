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


def get_diagnosis_service() -> DiagnosisService:
    """Dependency injection para DiagnosisService"""
    return DiagnosisService()


def get_semantic_search_service() -> SemanticSearchService:
    """Dependency injection para SemanticSearchService"""
    return SemanticSearchService()


def get_pubmed_service() -> PubMedService:
    """Dependency injection para PubMedService"""
    return PubMedService()


@diagnosis_router.post("/diagnosis", response_model=DiagnosisResponse)
async def create_diagnosis(
    request: DiagnosisRequest,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> DiagnosisResponse:
    """
    Crea un nuevo diagnóstico basado en síntomas y contexto médico
    """
    try:
        return await diagnosis_service.generate_diagnosis(request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al generar diagnóstico: {str(e)}")


@diagnosis_router.get("/diagnosis/history")
async def get_diagnosis_history(
    patient_id: str = None,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> List[Dict[str, Any]]:
    """
    Obtiene el historial de diagnósticos (opcionalmente filtrado por paciente)
    """
    try:
        return await diagnosis_service.get_diagnosis_history(patient_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al obtener historial: {str(e)}")


@diagnosis_router.post("/diagnosis/update-knowledge")
async def update_medical_knowledge(
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> Dict[str, Any]:
    """
    Actualiza la base de conocimiento médico
    """
    try:
        return await diagnosis_service.update_medical_knowledge()
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al actualizar conocimiento: {str(e)}")


@diagnosis_router.get("/health")
async def health_check() -> Dict[str, str]:
    """
    Endpoint de verificación de salud de la API
    """
    return {"status": "healthy", "service": "PySIO AI Diagnosis API"}


@diagnosis_router.post("/test/pubmed-search")
async def test_pubmed_search(
    query: str = "headache",
    max_results: int = 3,
    pubmed_service: PubMedService = Depends(get_pubmed_service)
) -> Dict[str, Any]:
    """
    Endpoint de prueba para búsqueda en PubMed
    """
    try:
        articles = await pubmed_service.search_and_fetch(query, max_results)
        return {
            "query": query,
            "articles_found": len(articles),
            "articles": articles[:2]  # Solo retornar los primeros 2 para no sobrecargar
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error en búsqueda PubMed: {str(e)}")


@diagnosis_router.get("/test/vector-store-stats")
async def test_vector_store_stats(
    semantic_search_service: SemanticSearchService = Depends(get_semantic_search_service)
) -> Dict[str, Any]:
    """
    Endpoint de prueba para estadísticas del vector store
    """
    try:
        stats = await semantic_search_service.get_vector_store_stats()
        return stats
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error obteniendo estadísticas: {str(e)}")


@diagnosis_router.post("/test/ingest-pubmed")
async def test_ingest_pubmed(
    query: str = "headache diagnosis",
    max_results: int = 3,
    semantic_search_service: SemanticSearchService = Depends(get_semantic_search_service)
) -> Dict[str, Any]:
    """
    Endpoint de prueba para ingestar artículos de PubMed
    """
    try:
        success = await semantic_search_service.ingest_pubmed_articles(query, max_results)
        return {
            "query": query,
            "ingestion_success": success,
            "message": "Artículos ingresados en vector store" if success else "Error en ingesta"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error en ingesta: {str(e)}")

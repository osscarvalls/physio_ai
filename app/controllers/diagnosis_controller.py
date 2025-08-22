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

logger = logging.getLogger(__name__)

# Router para los endpoints de diagnóstico
diagnosis_router = APIRouter(prefix="/api/v1", tags=["diagnosis"])


def get_diagnosis_service() -> DiagnosisService:
    """Dependency injection para DiagnosisService"""
    return DiagnosisService()


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

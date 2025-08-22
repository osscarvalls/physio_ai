"""
Controlador para endpoints de diagnóstico médico
"""

import logging
from fastapi import APIRouter, HTTPException, Depends
from fastapi.responses import HTMLResponse
from fastapi.requests import Request

from app.models.diagnosis import DiagnosisRequest, DiagnosisResponse
from app.services.diagnosis_service import DiagnosisService
from app.config.settings import settings

logger = logging.getLogger(__name__)

# Crear router para diagnóstico
diagnosis_router = APIRouter(prefix="/api/v1", tags=["diagnóstico"])


def get_diagnosis_service() -> DiagnosisService:
    """Dependency injection para el servicio de diagnóstico"""
    return DiagnosisService()


@diagnosis_router.get("/", response_class=HTMLResponse)
async def home(request: Request):
    """Página principal de la aplicación"""
    from fastapi.templating import Jinja2Templates
    
    templates = Jinja2Templates(directory=settings.TEMPLATES_DIR)
    return templates.TemplateResponse("index.html", {"request": request})


@diagnosis_router.post("/diagnosis", response_model=DiagnosisResponse)
async def create_diagnosis(
    request: DiagnosisRequest,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
) -> DiagnosisResponse:
    """
    Genera un diagnóstico médico basado en los síntomas proporcionados
    
    Args:
        request: Solicitud de diagnóstico con síntomas y datos del paciente
        diagnosis_service: Servicio de diagnóstico inyectado
        
    Returns:
        Respuesta con diagnóstico y recomendaciones
        
    Raises:
        HTTPException: Si hay error en el procesamiento
    """
    try:
        logger.info(f"Solicitud de diagnóstico recibida para síntomas: {request.symptoms[:100]}...")
        
        # Validar entrada
        if not request.symptoms or not request.symptoms.strip():
            raise HTTPException(
                status_code=400, 
                detail="Los síntomas son obligatorios"
            )
        
        # Generar diagnóstico
        response = await diagnosis_service.generate_diagnosis(request)
        
        logger.info("Diagnóstico generado exitosamente")
        return response
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error inesperado en endpoint de diagnóstico: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Error interno del servidor al generar diagnóstico"
        )


@diagnosis_router.get("/diagnosis/history")
async def get_diagnosis_history(
    patient_id: str = None,
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
):
    """
    Obtiene historial de diagnósticos previos
    
    Args:
        patient_id: ID del paciente (opcional)
        diagnosis_service: Servicio de diagnóstico inyectado
        
    Returns:
        Lista de diagnósticos previos
    """
    try:
        history = await diagnosis_service.get_diagnosis_history(patient_id)
        return {
            "status": "success",
            "data": history,
            "count": len(history)
        }
        
    except Exception as e:
        logger.error(f"Error obteniendo historial: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Error obteniendo historial de diagnósticos"
        )


@diagnosis_router.post("/diagnosis/update-knowledge")
async def update_medical_knowledge(
    diagnosis_service: DiagnosisService = Depends(get_diagnosis_service)
):
    """
    Actualiza la base de conocimiento médica
    
    Args:
        diagnosis_service: Servicio de diagnóstico inyectado
        
    Returns:
        Estado de la actualización
    """
    try:
        result = await diagnosis_service.update_medical_knowledge()
        return result
        
    except Exception as e:
        logger.error(f"Error actualizando conocimiento médico: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail="Error actualizando base de conocimiento médica"
        )


@diagnosis_router.get("/health")
async def health_check():
    """Endpoint de verificación de salud del servicio"""
    return {
        "status": "healthy",
        "service": "PySIO AI Diagnosis Service",
        "version": "1.0.0",
        "timestamp": "now"
    }

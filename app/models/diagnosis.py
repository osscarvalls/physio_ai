"""
Modelos de datos para el sistema de diagnóstico médico
"""

from datetime import datetime
from typing import List, Optional
from pydantic import BaseModel, Field

class DiagnosisRequest(BaseModel):
    """Modelo de solicitud para diagnóstico"""
    symptoms: str = Field(..., description="Descripción de los síntomas del paciente")
    patient_age: Optional[int] = Field(None, description="Edad del paciente")
    patient_gender: Optional[str] = Field(None, description="Género del paciente")
    medical_history: Optional[str] = Field(None, description="Historial médico relevante")


class DiagnosisResponse(BaseModel):
    """Modelo de respuesta de diagnóstico"""
    patient_situation: str = Field(..., description="Situación del paciente")
    diagnostic_suggestions: Optional[List[str]] = Field(None, description="Sugerencias de diagnóstico")
    confirmation_tests: Optional[List[str]] = Field(None, description="Pruebas de confirmación del diagnóstico")
    diagnosis: str = Field(..., description="Diagnóstico generado por la IA")
    confidence: Optional[float] = Field(None, description="Nivel de confianza del diagnóstico")
    evidence_quality: Optional[str] = Field(None, description="Valoración de la relevancia y calidad de la evidencia")
    missing_information: Optional[List[str]] = Field(None, description="Información que falta en la evidencia para generar un diagnóstico completo")
    recommendations: Optional[List[str]] = Field(None, description="Recomendaciones médicas")
    timestamp: datetime = Field(default_factory=datetime.now, description="Timestamp del diagnóstico")
    error: Optional[str] = Field(None, description="Mensaje de error si ocurre algún problema")
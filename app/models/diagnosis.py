"""
Modelos de datos para el sistema de diagnóstico médico
"""

from datetime import datetime
from typing import List, Optional
from pydantic import BaseModel, Field


class Symptom(BaseModel):
    """Modelo para representar un síntoma"""
    description: str = Field(..., description="Descripción del síntoma")
    severity: Optional[str] = Field(None, description="Severidad del síntoma")
    duration: Optional[str] = Field(None, description="Duración del síntoma")


class DiagnosisRequest(BaseModel):
    """Modelo de solicitud para diagnóstico"""
    symptoms: str = Field(..., description="Descripción de los síntomas del paciente")
    patient_age: Optional[int] = Field(None, description="Edad del paciente")
    patient_gender: Optional[str] = Field(None, description="Género del paciente")
    medical_history: Optional[str] = Field(None, description="Historial médico relevante")


class DiagnosisResponse(BaseModel):
    """Modelo de respuesta de diagnóstico"""
    diagnosis: str = Field(..., description="Diagnóstico generado por la IA")
    confidence: Optional[float] = Field(None, description="Nivel de confianza del diagnóstico")
    recommendations: Optional[List[str]] = Field(None, description="Recomendaciones médicas")
    timestamp: datetime = Field(default_factory=datetime.now, description="Timestamp del diagnóstico")


class MedicalEvidence(BaseModel):
    """Modelo para evidencia médica de PubMed"""
    pmid: str = Field(..., description="ID de PubMed")
    title: str = Field(..., description="Título del artículo")
    abstract: str = Field(..., description="Resumen del artículo")
    authors: Optional[List[str]] = Field(None, description="Autores del artículo")
    publication_date: Optional[str] = Field(None, description="Fecha de publicación")
    relevance_score: Optional[float] = Field(None, description="Puntuación de relevancia")


class SearchQuery(BaseModel):
    """Modelo para consultas de búsqueda"""
    query: str = Field(..., description="Consulta de búsqueda")
    symptoms: str = Field(..., description="Síntomas originales")
    generated_at: datetime = Field(default_factory=datetime.now, description="Timestamp de generación")

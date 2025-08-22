"""
Servicio principal de diagnóstico médico
"""

import logging
from typing import List, Dict, Any
from langchain.schema import Document

from app.services.llm_service import LLMService
from app.services.evidence_service import EvidenceService
from app.models.diagnosis import DiagnosisRequest, DiagnosisResponse, MedicalEvidence

logger = logging.getLogger(__name__)


class DiagnosisService:
    """Servicio principal para generar diagnósticos médicos"""
    
    def __init__(self):
        """Inicializa el servicio de diagnóstico"""
        self.llm_service = LLMService()
        self.evidence_service = EvidenceService()
    
    async def generate_diagnosis(self, request: DiagnosisRequest) -> DiagnosisResponse:
        """
        Genera un diagnóstico completo basado en la solicitud
        
        Args:
            request: Solicitud de diagnóstico con síntomas y datos del paciente
            
        Returns:
            Respuesta con diagnóstico y recomendaciones
        """
        try:
            logger.info(f"Generando diagnóstico para síntomas: {request.symptoms[:100]}...")
            
            # Paso 1: Generar consultas de búsqueda optimizadas
            search_queries = await self.llm_service.generate_search_queries(request.symptoms)
            logger.info(f"Consultas de búsqueda generadas: {len(search_queries)}")
            
            # Paso 2: Buscar evidencia médica relevante
            medical_evidence = await self._gather_medical_evidence(search_queries)
            logger.info(f"Evidencia médica encontrada: {len(medical_evidence)} documentos")
            
            # Paso 3: Generar diagnóstico usando LLM
            diagnosis_result = await self.llm_service.generate_diagnosis(
                symptoms=request.symptoms,
                medical_context=medical_evidence
            )
            
            # Paso 4: Construir respuesta
            response = DiagnosisResponse(
                diagnosis=diagnosis_result.get("diagnosis", "No se pudo generar diagnóstico"),
                confidence=diagnosis_result.get("confidence", 0.0),
                recommendations=self._generate_recommendations(diagnosis_result),
                timestamp=diagnosis_result.get("timestamp")
            )
            
            logger.info("Diagnóstico generado exitosamente")
            return response
            
        except Exception as e:
            logger.error(f"Error generando diagnóstico: {str(e)}")
            return DiagnosisResponse(
                diagnosis="No se pudo generar diagnóstico en este momento. Por favor, intente más tarde.",
                confidence=0.0,
                recommendations=["Consulte con un profesional médico"],
                error=str(e)
            )
    
    async def _gather_medical_evidence(self, search_queries: List[str]) -> List[Document]:
        """
        Recopila evidencia médica de múltiples fuentes
        
        Args:
            search_queries: Lista de consultas de búsqueda
            
        Returns:
            Lista de documentos médicos relevantes
        """
        all_evidence = []
        
        for query in search_queries[:3]:  # Limitar a 3 consultas para evitar sobrecarga
            try:
                # Buscar en PubMed
                pmid_list = await self.evidence_service.search_pubmed(query, max_results=5)
                
                if pmid_list:
                    # Obtener detalles de los artículos
                    papers = await self.evidence_service.fetch_papers(pmid_list)
                    
                    # Procesar y almacenar localmente
                    if papers:
                        await self.evidence_service.store_papers_locally(papers)
                        
                        # Convertir a documentos de LangChain
                        for paper in papers:
                            try:
                                pmid = paper['MedlineCitation']['PMID']
                                title = paper['MedlineCitation']['Article']['ArticleTitle']
                                
                                # Obtener abstract
                                abstract = ""
                                if 'Abstract' in paper['MedlineCitation']['Article']:
                                    abstract_text = paper['MedlineCitation']['Article']['Abstract'].get('AbstractText', [])
                                    if isinstance(abstract_text, list):
                                        abstract = " ".join(abstract_text)
                                    else:
                                        abstract = str(abstract_text)
                                
                                # Crear documento de LangChain
                                doc = Document(
                                    page_content=f"Title: {title}\nAbstract: {abstract}",
                                    metadata={
                                        "pmid": pmid,
                                        "title": title,
                                        "source": "pubmed",
                                        "query": query
                                    }
                                )
                                all_evidence.append(doc)
                                
                            except Exception as e:
                                logger.error(f"Error procesando artículo: {str(e)}")
                                continue
                
            except Exception as e:
                logger.error(f"Error procesando consulta '{query}': {str(e)}")
                continue
        
        # Si no hay evidencia de PubMed, usar documentos locales
        if not all_evidence:
            logger.info("No se encontró evidencia de PubMed, usando documentos locales")
            all_evidence = await self.evidence_service.retrieve_relevant_evidence(
                symptoms=" ".join(search_queries),
                max_results=5
            )
        
        return all_evidence
    
    def _generate_recommendations(self, diagnosis_result: Dict[str, Any]) -> List[str]:
        """
        Genera recomendaciones médicas basadas en el diagnóstico
        
        Args:
            diagnosis_result: Resultado del diagnóstico del LLM
            
        Returns:
            Lista de recomendaciones
        """
        base_recommendations = [
            "Consulte con un profesional médico para confirmar el diagnóstico",
            "Mantenga un registro de sus síntomas y su evolución",
            "No se automedique sin supervisión médica"
        ]
        
        # Agregar recomendaciones específicas si hay confianza alta
        if diagnosis_result.get("confidence", 0.0) > 0.7:
            base_recommendations.append(
                "El diagnóstico tiene alta confianza, pero aún requiere confirmación médica"
            )
        
        # Agregar recomendaciones basadas en el contexto usado
        context_used = diagnosis_result.get("context_used", 0)
        if context_used > 0:
            base_recommendations.append(
                f"El diagnóstico se basó en {context_used} fuentes médicas relevantes"
            )
        
        return base_recommendations
    
    async def update_medical_knowledge(self) -> Dict[str, Any]:
        """
        Actualiza la base de conocimiento médica
        
        Returns:
            Dict con información sobre la actualización
        """
        try:
            logger.info("Iniciando actualización de base de conocimiento médica")
            
            # Actualizar vector store
            await self.evidence_service.update_vector_store()
            
            return {
                "status": "success",
                "message": "Base de conocimiento médica actualizada",
                "timestamp": "now"
            }
            
        except Exception as e:
            logger.error(f"Error actualizando base de conocimiento: {str(e)}")
            return {
                "status": "error",
                "message": f"Error actualizando base de conocimiento: {str(e)}",
                "timestamp": "now"
            }
    
    async def get_diagnosis_history(self, patient_id: str = None) -> List[Dict[str, Any]]:
        """
        Obtiene historial de diagnósticos (placeholder para futura implementación)
        
        Args:
            patient_id: ID del paciente (opcional)
            
        Returns:
            Lista de diagnósticos previos
        """
        # TODO: Implementar persistencia de diagnósticos
        logger.info("Historial de diagnósticos solicitado (no implementado)")
        return []

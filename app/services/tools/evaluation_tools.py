"""
Herramientas de evaluación para el servicio de diagnóstico
"""

import logging
import json
from typing import List, Dict, Any
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI

logger = logging.getLogger(__name__)


class EvaluationTools:
    """Herramientas para evaluar la relevancia de la evidencia médica"""
    
    def __init__(self, llm: ChatOpenAI):
        self.llm = llm
    
    async def evaluate_relevance(self, symptoms: str, medical_evidence: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Evalúa la relevancia de la evidencia encontrada"""
        try:
            # Si no hay evidencia, evaluar como no relevante
            if not medical_evidence:
                return {
                    "relevance_score": 0,
                    "is_relevant": False,
                    "reason": "No se encontró evidencia médica",
                    "missing_aspects": ["evidencia médica"]
                }
            
            # Formatear contexto médico
            context = self._format_medical_context(medical_evidence)
            
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta experto que evalúa la relevancia de la evidencia médica.
            
            Síntomas del paciente: {symptoms}
            Evidencia médica disponible: {context}
            
            Evalúa si la evidencia es RELEVANTE para generar un diagnóstico:
            
            CRITERIOS:
            - La evidencia debe estar directamente relacionada con los síntomas
            - Debe incluir información sobre diagnóstico, evaluación o tratamiento
            - Debe ser específica para fisioterapia
            - Debe ser de calidad suficiente
            
            Responde SOLO con este JSON:
            {{
                "relevance_score": <0-10>,
                "is_relevant": <true/false>,
                "reason": "explicación de la evaluación",
                "missing_aspects": ["aspecto1", "aspecto2"]
            }}
            """)
            
            chain = prompt | self.llm
            result = await chain.ainvoke({
                "symptoms": symptoms,
                "context": context
            })
            
            # Parsear resultado JSON
            try:
                # Limpiar el contenido del LLM (eliminar markdown si existe)
                content = result.content.strip()
                if content.startswith('```json'):
                    content = content[7:]  # Eliminar ```json
                if content.endswith('```'):
                    content = content[:-3]  # Eliminar ```
                content = content.strip()
                
                evaluation = json.loads(content)
                logger.info(f"Evaluación de relevancia: {evaluation}")
                return evaluation
            except json.JSONDecodeError:
                logger.error("Error parseando evaluación de relevancia")
                return {
                    "relevance_score": 3,
                    "is_relevant": False,
                    "reason": "Error en evaluación automática",
                    "missing_aspects": ["evaluación manual requerida"]
                }
            
        except Exception as e:
            logger.error(f"Error evaluando relevancia: {str(e)}")
            return {
                "relevance_score": 0,
                "is_relevant": False,
                "reason": f"Error en evaluación: {str(e)}",
                "missing_aspects": ["evaluación manual requerida"]
            }
    
    def evaluate_relevance_router(self, evaluation: Dict[str, Any], iterations: int, max_iterations: int) -> str:
        """Determina el siguiente paso basado en la evaluación de relevancia"""
        try:
            # Si ya alcanzamos el máximo de iteraciones
            if iterations >= max_iterations:
                logger.info(f"Alcanzado máximo de iteraciones ({max_iterations})")
                return "no_evidence"
            
            # Si la evidencia es relevante, proceder con diagnóstico
            if evaluation.get("is_relevant", False):
                logger.info("Evidencia relevante encontrada, procediendo con diagnóstico")
                return "relevant"
            
            # Si no es relevante y podemos iterar más, expandir bibliografía
            logger.info(f"Evidencia no relevante, expandiendo bibliografía (iteración {iterations + 1})")
            return "expand"
            
        except Exception as e:
            logger.error(f"Error en lógica condicional: {str(e)}")
            return "expand"  # Por defecto, expandir
    
    def _format_medical_context(self, evidence: List[Dict[str, Any]]) -> str:
        """Formatea la evidencia médica para el prompt"""
        if not evidence:
            return "No se encontró evidencia médica"
        
        context_parts = []
        for i, doc in enumerate(evidence[:5], 1):  # Limitar a 5 documentos
            title = doc.get('title', 'Sin título')
            content = doc.get('content', 'Sin contenido')
            context_parts.append(f"Documento {i}: {title}\n{content[:200]}...")
        
        return "\n\n".join(context_parts)

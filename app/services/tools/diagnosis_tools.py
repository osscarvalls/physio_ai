"""
Herramientas de diagnóstico para el servicio de diagnóstico
"""

import logging
import json
from typing import List, Dict, Any
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI

logger = logging.getLogger(__name__)


class DiagnosisTools:
    """Herramientas para generar diagnósticos y evaluaciones del paciente"""
    
    def __init__(self, llm: ChatOpenAI):
        self.llm = llm
    
    async def generate_diagnosis(self, symptoms: str, age: int, gender: str, 
                                medical_evidence: List[Dict[str, Any]], 
                                evidence_evaluation: Dict[str, Any]) -> Dict[str, Any]:
        """Genera la evaluación completa del paciente para el fisioterapeuta"""
        try:
            # Formatear contexto médico
            context = self._format_medical_context(medical_evidence)
            
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta experto. Genera una evaluación del paciente en formato JSON.
            
            PACIENTE: {age} años, {gender}, síntomas: {symptoms}
            EVIDENCIA MÉDICA: {context}
            EVALUACIÓN DE RELEVANCIA: {evaluation}
            
            Responde SOLO con este JSON exacto:
            {{
                "patient_situation": "Descripción de la situación del paciente",
                "diagnosis_summary": "Resumen del diagnóstico basado en la evidencia",
                "confidence": "Confianza en el diagnóstico, entre 0 y 1", 
                "confirmation_tests": ["Prueba 1", "Prueba 2"], (pruebas de confirmación del diagnóstico)
                "recommendations": ["Recomendación 1", "Recomendación 2"], (recomendaciones para el paciente si encaja con el diagnóstico)
                "evidence_quality": "Valoración de la relevancia y calidad de la evidencia", (valoración de la relevancia y calidad de la evidencia)
                "missing_information": ["Información faltante 1", "Información faltante 2"] (información que falta en la evidencia para generar un diagnóstico completo)
            }}

            NOTA IMPORTANTE: TODO lo que escribas en el JSON debe estar basado en la evidencia y los síntomas del paciente. 
            No inventes información. Si hay algo que no esté claro, escribe "No se puede determinar" o "No se puede generar un diagnóstico completo".
            Incluye en el campo Missing Information la información que falta para generar un diagnóstico completo. 
            Esta información volverá a un evaluador que decidirá cómo ampliar la evidencia para generar un diagnóstico completo.
            """)
            
            chain = prompt | self.llm
            result = await chain.ainvoke({
                "age": age,
                "gender": gender,
                "symptoms": symptoms,
                "context": context,
                "evaluation": evidence_evaluation
            })
            
            logger.info(f"Respuesta del LLM: {result.content}")
            
            # Parsear resultado JSON
            try:
                # Limpiar el contenido del LLM (eliminar markdown si existe)
                content = result.content.strip()
                if content.startswith('```json'):
                    content = content[7:] 
                if content.endswith('```'):
                    content = content[:-3]
                content = content.strip()
                
                patient_evaluation = json.loads(content)
                logger.info("Evaluación del paciente generada exitosamente")
                return patient_evaluation

            except json.JSONDecodeError as e:
                logger.error(f"Error parseando evaluación del paciente: {e}")
                logger.error(f"Contenido del LLM: {result.content}")
                return {
                    "diagnosis_summary": "No se pudo generar evaluación completa",
                    "confidence": 0.5,
                    "patient_situation": "Paciente con síntomas que requieren evaluación",
                    "diagnostic_suggestions": ["Evaluación fisioterapéutica completa"],
                    "confirmation_tests": ["Exploración física", "Historia clínica"],
                    "recommendations": ["Consulta con fisioterapeuta"],
                    "evidence_quality": "baja",
                    "missing_information": ["Evidencia insuficiente"],
                    "evidence_based": False
                }
            
        except Exception as e:
            logger.error(f"Error generando evaluación del paciente: {str(e)}")
            return {
                "diagnosis_summary": "Error generando evaluación",
                "confidence": 0.0,
                "patient_situation": "Error en el sistema",
                "diagnostic_suggestions": ["Evaluación manual requerida"],
                "confirmation_tests": ["Revisión del sistema"],
                "recommendations": ["Contactar soporte técnico"],
                "evidence_quality": "error",
                "missing_information": ["Funcionamiento del sistema"],
                "evidence_based": False
            }
    
    def generate_no_evidence_response(self) -> Dict[str, Any]:
        """Genera una respuesta cuando no se encuentra evidencia médica"""
        try:
            return {
                "diagnosis_summary": "No se encontró evidencia médica suficiente para generar un diagnóstico",
                "confidence": 0.0,
                "patient_situation": "Paciente con síntomas que requieren evaluación",
                "diagnostic_suggestions": ["Evaluación fisioterapéutica completa"],
                "confirmation_tests": ["Exploración física", "Historia clínica"],
                "recommendations": ["Consulta con fisioterapeuta"],
                "evidence_quality": "baja",
                "missing_information": ["Evidencia insuficiente"],
                "evidence_based": False
            }
        except Exception as e:
            logger.error(f"Error generando respuesta cuando no se encuentra evidencia médica: {str(e)}")
            return {
                "diagnosis_summary": "Error en el sistema",
                "confidence": 0.0,
                "patient_situation": "Error en el sistema",
                "diagnostic_suggestions": ["Evaluación manual requerida"],
                "confirmation_tests": ["Revisión del sistema"],
                "recommendations": ["Contactar soporte técnico"],
                "evidence_quality": "error",
                "missing_information": ["Funcionamiento del sistema"],
                "evidence_based": False
            }
    
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

"""
Servicio principal de diagnóstico médico con LangGraph
"""

import logging
from typing import List, Dict, Any, TypedDict, Annotated
from datetime import datetime
from langgraph.graph import StateGraph, START, END
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

from app.services.semantic_search_service import SemanticSearchService
from app.models.diagnosis import DiagnosisRequest, DiagnosisResponse, MedicalEvidence
from app.config.settings import settings

logger = logging.getLogger(__name__)


# Estados del grafo
class DiagnosisState(TypedDict):
    """Estado del flujo de diagnóstico"""
    symptoms: str
    patient_age: int
    patient_gender: str
    search_queries: List[str]
    medical_evidence: List[Dict[str, Any]]
    evidence_evaluation: Dict[str, Any]
    search_iterations: int
    max_iterations: int
    patient_evaluation: Dict[str, Any]
    error: str
    timestamp: str


class DiagnosisService:
    """Servicio principal para generar diagnósticos médicos con LangGraph"""
    
    def __init__(self):
        """Inicializa el servicio de diagnóstico"""
        self.semantic_search_service = SemanticSearchService()
        self.llm = ChatOpenAI(
            model='gpt-4o-mini',
            temperature=0.1,
            api_key=settings.OPENAI_API_KEY
        )
        self._setup_graph()
    
    def _setup_graph(self):
        """Configura el grafo de LangGraph"""
        # Crear el grafo
        workflow = StateGraph(DiagnosisState)
        
        # Añadir nodos
        workflow.add_node("generate_queries", self._generate_search_queries)
        workflow.add_node("generate_patient_evaluation", self._generate_patient_evaluation)
        
        # Definir el flujo simple (solo 2 nodos por ahora)
        workflow.add_edge(START, "generate_queries")
        workflow.add_edge("generate_queries", "generate_patient_evaluation")
        workflow.add_edge("generate_patient_evaluation", END)
        
        # Compilar el grafo
        self.graph = workflow.compile()
        
        logger.info("Grafo de diagnóstico configurado con LangGraph (versión mínima)")
    
    async def generate_diagnosis(self, request: DiagnosisRequest) -> DiagnosisResponse:
        """
        Genera una evaluación completa del paciente basada en la solicitud
        
        Args:
            request: Solicitud de diagnóstico con síntomas y datos del paciente
            
        Returns:
            Respuesta con evaluación del paciente
        """
        try:
            logger.info(f"Generando evaluación para síntomas: {request.symptoms[:100]}...")
            
            # Estado inicial
            initial_state = DiagnosisState(
                symptoms=request.symptoms,
                patient_age=request.patient_age,
                patient_gender=request.patient_gender,
                search_queries=[],
                medical_evidence=[],
                evidence_evaluation={},
                search_iterations=0,
                max_iterations=3,
                patient_evaluation={},
                error="",
                timestamp=datetime.now().isoformat()
            )
            
            logger.info(f"Estado inicial creado: {initial_state}")
            
            # Ejecutar el grafo
            final_state = await self.graph.ainvoke(initial_state)
            
            logger.info(f"Estado final del grafo: {final_state}")
            
            # Construir respuesta
            patient_eval = final_state.get("patient_evaluation", {})
            response = DiagnosisResponse(
                diagnosis=patient_eval.get("diagnosis_summary", "No se pudo generar evaluación"),
                confidence=patient_eval.get("confidence", 0.0),
                recommendations=patient_eval.get("recommendations", []),
                timestamp=final_state.get("timestamp")
            )
            
            logger.info("Evaluación del paciente generada exitosamente")
            return response
            
        except Exception as e:
            logger.error(f"Error generando evaluación: {str(e)}")
            return DiagnosisResponse(
                diagnosis="No se pudo generar evaluación en este momento. Por favor, intente más tarde.",
                confidence=0.0,
                recommendations=["Consulte con un profesional médico"],
                error=str(e)
            )
    
    async def _generate_search_queries(self, state: DiagnosisState) -> DiagnosisState:
        """Genera consultas de búsqueda optimizadas para fisioterapia"""
        try:
            symptoms = state["symptoms"]
            age = state["patient_age"]
            gender = state["patient_gender"]
            
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta experto que necesita buscar evidencia médica para evaluar un paciente.
            
            Paciente: {age} años, {gender}
            Síntomas: {symptoms}
            
            Genera 3-5 consultas de búsqueda optimizadas para PubMed que se enfoquen en:
            1. Diagnóstico diferencial de fisioterapia
            2. Evaluación y pruebas específicas
            3. Tratamiento fisioterapéutico
            4. Evidencia clínica relevante
            
            Las consultas deben ser específicas y precisas, priorizando la calidad sobre la cantidad.
            Genera solo las consultas, una por línea, sin numeración ni formato adicional.
            """)
            
            chain = prompt | self.llm
            result = await chain.ainvoke({
                "symptoms": symptoms,
                "age": age,
                "gender": gender
            })
            
            # Parsear resultado
            queries = [
                query.strip() 
                for query in result.content.split('\n') 
                if query.strip()
            ]
            
            state["search_queries"] = queries[:5]  # Limitar a 5 consultas
            logger.info(f"Consultas de búsqueda generadas: {len(queries)}")
            
            return state
            
        except Exception as e:
            logger.error(f"Error generando consultas: {str(e)}")
            state["error"] = f"Error generando consultas: {str(e)}"
            return state
    
    async def _gather_medical_evidence(self, state: DiagnosisState) -> DiagnosisState:
        """Recopila evidencia médica relevante"""
        try:
            search_queries = state["search_queries"]
            all_evidence = []
            
            # Incrementar contador de iteraciones
            state["search_iterations"] += 1
            logger.info(f"Iteración de búsqueda {state['search_iterations']}/{state['max_iterations']}")
            
            for query in search_queries[:3]:  # Limitar a 3 consultas
                try:
                    # Buscar evidencia médica
                    evidence = await self.semantic_search_service.search_medical_evidence(query, max_results=3)
                    
                    if evidence:
                        all_evidence.extend(evidence)
                        logger.info(f"Encontrada evidencia para consulta '{query}': {len(evidence)} documentos")
                    else:
                        # Si no hay evidencia, intentar ingestar artículos de PubMed
                        logger.info(f"No hay evidencia para '{query}', intentando ingestar de PubMed")
                        success = await self.semantic_search_service.ingest_pubmed_articles(query, max_results=3)
                        if success:
                            # Buscar nuevamente
                            evidence = await self.semantic_search_service.search_medical_evidence(query, max_results=3)
                            if evidence:
                                all_evidence.extend(evidence)
                
                except Exception as e:
                    logger.error(f"Error procesando consulta '{query}': {str(e)}")
                    continue
            
            # Si no hay evidencia, usar búsqueda local
            if not all_evidence:
                logger.info("No se encontró evidencia, usando búsqueda local")
                all_evidence = await self.semantic_search_service.search_medical_evidence(
                    " ".join(search_queries),
                    max_results=5
                )
            
            state["medical_evidence"] = all_evidence
            logger.info(f"Evidencia médica encontrada: {len(all_evidence)} documentos")
            
            return state
            
        except Exception as e:
            logger.error(f"Error recopilando evidencia: {str(e)}")
            state["error"] = f"Error recopilando evidencia: {str(e)}"
            return state
    
    async def _evaluate_evidence(self, state: DiagnosisState) -> DiagnosisState:
        """Evalúa la calidad y cantidad de evidencia recopilada"""
        try:
            medical_evidence = state["medical_evidence"]
            symptoms = state["symptoms"]
            
            # Formatear contexto médico
            context = self._format_medical_context(medical_evidence)
            
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta experto que evalúa la calidad de la evidencia médica recopilada.
            
            Síntomas del paciente: {symptoms}
            Evidencia médica disponible: {context}
            
            Evalúa la evidencia en términos de:
            1. Relevancia para el caso (0-10)
            2. Calidad de la evidencia (0-10)
            3. Cobertura de aspectos clave (0-10)
            4. Suficiencia para generar una evaluación completa (SÍ/NO)
            
            Responde en formato JSON:
            {{
                "relevance_score": <puntuación>,
                "quality_score": <puntuación>,
                "coverage_score": <puntuación>,
                "is_sufficient": <SÍ/NO>,
                "missing_aspects": ["aspecto1", "aspecto2"],
                "recommendations": ["recomendación1", "recomendación2"]
            }}
            """)
            
            chain = prompt | self.llm
            result = await chain.ainvoke({
                "symptoms": symptoms,
                "context": context
            })
            
            # Parsear resultado JSON
            try:
                import json
                evaluation = json.loads(result.content)
                state["evidence_evaluation"] = evaluation
                logger.info(f"Evaluación de evidencia: {evaluation}")
            except json.JSONDecodeError:
                logger.error("Error parseando evaluación de evidencia")
                state["evidence_evaluation"] = {
                    "relevance_score": 5,
                    "quality_score": 5,
                    "coverage_score": 5,
                    "is_sufficient": "NO",
                    "missing_aspects": ["evidencia insuficiente"],
                    "recommendations": ["ampliar búsqueda"]
                }
            
            return state
            
        except Exception as e:
            logger.error(f"Error evaluando evidencia: {str(e)}")
            state["error"] = f"Error evaluando evidencia: {str(e)}"
            return state
    
    def _should_expand_evidence(self, state: DiagnosisState) -> str:
        """Determina si se debe expandir la evidencia o proceder con la evaluación"""
        try:
            evaluation = state.get("evidence_evaluation", {})
            iterations = state.get("search_iterations", 0)
            max_iterations = state.get("max_iterations", 3)
            
            # Si ya alcanzamos el máximo de iteraciones, proceder
            if iterations >= max_iterations:
                logger.info(f"Alcanzado máximo de iteraciones ({max_iterations}), procediendo con evaluación")
                return "sufficient"
            
            # Si la evidencia es suficiente, proceder
            is_sufficient = evaluation.get("is_sufficient", "NO")
            if is_sufficient == "SÍ":
                logger.info("Evidencia suficiente encontrada, procediendo con evaluación")
                return "sufficient"
            
            # Si la evidencia no es suficiente y podemos iterar más, expandir
            logger.info(f"Evidencia insuficiente, expandiendo búsqueda (iteración {iterations + 1})")
            return "expand"
            
        except Exception as e:
            logger.error(f"Error en lógica condicional: {str(e)}")
            return "sufficient"  # Por defecto, proceder
    
    async def _expand_evidence(self, state: DiagnosisState) -> DiagnosisState:
        """Expande la evidencia buscando en PubMed y actualizando la base de datos vectorial"""
        try:
            symptoms = state["symptoms"]
            evaluation = state.get("evidence_evaluation", {})
            missing_aspects = evaluation.get("missing_aspects", [])
            
            logger.info(f"Expandiendo evidencia para aspectos faltantes: {missing_aspects}")
            
            # Generar consultas específicas para aspectos faltantes
            if missing_aspects:
                # Crear consultas específicas para aspectos faltantes
                specific_queries = []
                for aspect in missing_aspects[:2]:  # Limitar a 2 aspectos
                    query = f"{aspect} {symptoms} fisioterapia"
                    specific_queries.append(query)
                
                # Añadir a las consultas existentes
                state["search_queries"].extend(specific_queries)
                logger.info(f"Consultas adicionales generadas: {specific_queries}")
            
            # Buscar e ingestar nuevos artículos de PubMed
            for query in state["search_queries"][-2:]:  # Solo las 2 últimas consultas
                try:
                    success = await self.semantic_search_service.ingest_pubmed_articles(query, max_results=2)
                    if success:
                        logger.info(f"Evidencia expandida para consulta: {query}")
                except Exception as e:
                    logger.error(f"Error expandiendo evidencia para '{query}': {str(e)}")
                    continue
            
            logger.info("Expansión de evidencia completada")
            return state
            
        except Exception as e:
            logger.error(f"Error expandiendo evidencia: {str(e)}")
            state["error"] = f"Error expandiendo evidencia: {str(e)}"
            return state
    
    def _format_medical_context(self, medical_evidence: List[Dict[str, Any]]) -> str:
        """Formatea el contexto médico para el prompt"""
        if not medical_evidence:
            return "No hay evidencia médica disponible."
        
        context_parts = []
        for i, evidence in enumerate(medical_evidence[:5], 1):  # Limitar a 5 documentos
            content = evidence.get('content', '')[:500]
            context_parts.append(f"Documento {i}: {content}...")
        
        return "\n\n".join(context_parts)
    
    async def _generate_patient_evaluation(self, state: DiagnosisState) -> DiagnosisState:
        """Genera la evaluación completa del paciente para el fisioterapeuta"""
        try:
            symptoms = state["symptoms"]
            patient_age = state["patient_age"]
            patient_gender = state["patient_gender"]
            medical_evidence = state["medical_evidence"]
            evidence_evaluation = state.get("evidence_evaluation", {})
            
            # Formatear contexto médico
            context = self._format_medical_context(medical_evidence)
            
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta experto. Genera una evaluación del paciente en formato JSON.
            
            PACIENTE: {age} años, {gender}, síntomas: {symptoms}
            
            Responde SOLO con este JSON exacto:
            {{
                "diagnosis_summary": "Resumen del diagnóstico",
                "confidence": 0.8,
                "patient_situation": "Descripción de la situación",
                "diagnostic_suggestions": ["sugerencia1", "sugerencia2"],
                "confirmation_tests": ["prueba1", "prueba2"],
                "recommendations": ["recomendación1", "recomendación2"],
                "evidence_quality": "alta",
                "missing_information": ["info1"]
            }}
            """)
            
            chain = prompt | self.llm
            result = await chain.ainvoke({
                "age": patient_age,
                "gender": patient_gender,
                "symptoms": symptoms
            })
            
            logger.info(f"Respuesta del LLM: {result.content}")
            
            # Parsear resultado JSON
            try:
                import json
                patient_evaluation = json.loads(result.content)
                state["patient_evaluation"] = patient_evaluation
                logger.info("Evaluación del paciente generada exitosamente")
            except json.JSONDecodeError as e:
                logger.error(f"Error parseando evaluación del paciente: {e}")
                logger.error(f"Contenido del LLM: {result.content}")
                state["patient_evaluation"] = {
                    "diagnosis_summary": "No se pudo generar evaluación completa",
                    "confidence": 0.5,
                    "patient_situation": "Paciente con síntomas que requieren evaluación",
                    "diagnostic_suggestions": ["Evaluación fisioterapéutica completa"],
                    "confirmation_tests": ["Exploración física", "Historia clínica"],
                    "recommendations": ["Consulta con fisioterapeuta"],
                    "evidence_quality": "baja",
                    "missing_information": ["Evidencia insuficiente"]
                }
            
            return state
            
        except Exception as e:
            logger.error(f"Error generando evaluación del paciente: {str(e)}")
            state["error"] = f"Error generando evaluación del paciente: {str(e)}"
            return state
    
    async def update_medical_knowledge(self) -> Dict[str, Any]:
        """
        Actualiza la base de conocimiento médico
        
        Returns:
            Dict con información sobre la actualización
        """
        try:
            logger.info("Iniciando actualización de base de conocimiento médica")
            
            # Actualizar vector store con nuevos artículos
            success = await self.semantic_search_service.ingest_pubmed_articles(
                "medical diagnosis", 
                max_results=10
            )
            
            return {
                "status": "success" if success else "error",
                "message": "Base de conocimiento médica actualizada" if success else "Error en actualización",
                "timestamp": datetime.now().isoformat()
            }
            
        except Exception as e:
            logger.error(f"Error actualizando base de conocimiento: {str(e)}")
            return {
                "status": "error",
                "message": f"Error actualizando base de conocimiento: {str(e)}",
                "timestamp": datetime.now().isoformat()
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

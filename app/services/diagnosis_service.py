"""
Servicio principal de diagnóstico médico con LangGraph
"""

import logging
from typing import List, Dict, Any, TypedDict, Annotated
from datetime import datetime
from langgraph.graph import StateGraph, START, END
from langchain_openai import ChatOpenAI


from app.services.semantic_search_service import SemanticSearchService
from app.services.pubmed_service import PubMedService
from app.services.tools import SearchTools, EvaluationTools, DiagnosisTools
from app.models.diagnosis import DiagnosisRequest, DiagnosisResponse
from app.config.settings import settings

logger = logging.getLogger(__name__)


# Estados del grafo
class DiagnosisState(TypedDict):
    # Datos del paciente
    symptoms: str
    patient_age: int
    patient_gender: str

    # Datos de la búsqueda
    search_queries: List[str]

    # Datos de la evidencia
    medical_evidence: List[Dict[str, Any]]
    evidence_evaluation: str

    # Datos de la búsqueda en PubMed
    pubmed_queries: List[str]

    # Control de flujo de búsqueda
    search_iterations: int
    max_search_iterations: int

    # Control de flujo de evaluación
    evaluation_iterations: int
    max_evaluation_iterations: int

    # Datos de la evaluación del paciente
    patient_evaluation: Dict[str, Any]

    # Control de errores
    error: str
    timestamp: str
    
class DiagnosisService:
    """Servicio principal para generar diagnósticos médicos con LangGraph"""
    
    def __init__(self):
        """Inicializa el servicio de diagnóstico"""
        self.semantic_search_service = SemanticSearchService()
        self.pubmed_service = PubMedService()
        self.llm = ChatOpenAI(
            model='gpt-5-nano',
            temperature=0.1,
            api_key=settings.OPENAI_API_KEY
        )
        
        # Inicializar herramientas
        self.search_tools = SearchTools(self.llm)
        self.evaluation_tools = EvaluationTools(self.llm)
        self.diagnosis_tools = DiagnosisTools(self.llm)
        
        self._setup_graph()
    
    def _setup_graph(self):
        """Configura el grafo de LangGraph con flujo iterativo y condicional"""
        # Crear el grafo
        workflow = StateGraph(DiagnosisState)
        
        # Añadir nodos
        workflow.add_node("generate_semantic_queries", self._generate_semantic_queries)
        workflow.add_node("search_semantic", self._search_semantic)
        workflow.add_node("evaluate_relevance", self._evaluate_relevance)
        workflow.add_node("generate_pubmed_queries", self._generate_pubmed_queries)
        workflow.add_node("search_pubmed", self._search_pubmed)
        workflow.add_node("generate_diagnosis", self._generate_diagnosis)
        workflow.add_node("no_evidence_found", self._no_evidence_found)
        
        # Definir el flujo iterativo y condicional
        workflow.add_edge(START, "generate_semantic_queries")
        workflow.add_edge("generate_semantic_queries", "search_semantic")
        workflow.add_edge("search_semantic", "evaluate_relevance")
        
        # Rama condicional: ¿Es relevante la evidencia?
        workflow.add_conditional_edges(
            "evaluate_relevance",
            self._evaluate_relevance_router,
            {
                "relevant": "generate_diagnosis",
                "expand": "generate_pubmed_queries",
                "no_evidence": "no_evidence_found"
            }
        )
        
        # Rama de expansión: volver a buscar después de expandir
        workflow.add_edge("generate_pubmed_queries", "search_pubmed")
        workflow.add_edge("search_pubmed", "search_semantic")
        workflow.add_edge("search_semantic", "evaluate_relevance")
        
        # Nodos finales
        workflow.add_edge("generate_diagnosis", END)
        workflow.add_edge("no_evidence_found", END)
        
        # Compilar el grafo
        self.graph = workflow.compile()
        
        logger.info("Grafo de diagnóstico configurado con flujo iterativo y condicional")
        
    async def evaluate_patient(self, request: DiagnosisRequest) -> DiagnosisResponse:
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
                max_search_iterations=3,
                evaluation_iterations=0,
                max_evaluation_iterations=2,
                patient_evaluation={},
                error="",
                timestamp=datetime.now().isoformat()
            )
            
            logger.info(f"Estado inicial creado: {initial_state}")
            
            # Ejecutar el grafo
            final_state = await self.graph.ainvoke(initial_state)
            
            logger.info(f"Estado final del grafo: {final_state}")
            
            # Construir respuesta según el nodo final
            if "patient_evaluation" in final_state and final_state["patient_evaluation"]:
                patient_eval = final_state["patient_evaluation"]
                response = DiagnosisResponse(
                    diagnosis=patient_eval.get("diagnosis_summary", "No se pudo generar evaluación"),
                    confidence=patient_eval.get("confidence", 0.0),
                    recommendations=patient_eval.get("recommendations", []),
                    timestamp=final_state.get("timestamp")
                )
                logger.info("Evaluación del paciente generada exitosamente")
            else:
                # Caso de no evidencia encontrada
                response = DiagnosisResponse(
                    diagnosis="No se encontró evidencia médica suficiente para generar un diagnóstico",
                    confidence=0.0,
                    recommendations=["Consulte con un profesional médico para evaluación directa"],
                    timestamp=final_state.get("timestamp")
                )
                logger.info("No se encontró evidencia médica suficiente")
            
            return response
            
        except Exception as e:
            logger.error(f"Error generando evaluación: {str(e)}")
            return DiagnosisResponse(
                diagnosis="No se pudo generar evaluación en este momento. Por favor, intente más tarde.",
                confidence=0.0,
                recommendations=["Consulte con un profesional médico"],
                error=str(e)
            )
    


    async def _generate_semantic_queries(self, state: DiagnosisState) -> DiagnosisState:
        """Genera consultas de búsqueda optimizadas para fisioterapia"""
        try:
            symptoms = state["symptoms"]
            age = state["patient_age"]
            gender = state["patient_gender"]
            
            queries = await self.search_tools.generate_semantic_queries(symptoms, age, gender)
            state["search_queries"] = queries
            
            if not queries:
                state["error"] = "Error generando consultas de búsqueda"
            
            return state
            
        except Exception as e:
            logger.error(f"Error generando consultas: {str(e)}")
            state["error"] = f"Error generando consultas: {str(e)}"
            return state
    
    async def _search_semantic(self, state: DiagnosisState) -> DiagnosisState:
        """Realiza búsqueda semántica en Qdrant"""
        try:
            search_queries = state["search_queries"]
            symptoms = state["symptoms"]
            all_evidence = []
            
            # Incrementar contador de iteraciones
            state["search_iterations"] += 1
            logger.info(f"Iteración de búsqueda {state['search_iterations']}/{state['max_search_iterations']}")
            
            # Validar que tengamos consultas de búsqueda
            if not search_queries:
                logger.warning("No hay consultas de búsqueda disponibles")
                state["medical_evidence"] = []
                return state
            
            # Buscar evidencia semántica para cada consulta
            for query in search_queries[:3]:  # Limitar a 3 consultas por iteración
                try:
                    evidence = await self.semantic_search_service.search_medical_evidence(query, max_results=3)
                    if evidence:
                        # Evitar duplicados por PMID antes de agregar a all_evidence
                        pmids_existentes = {doc.get('metadata', {}).get('pmid') for doc in all_evidence if doc.get('metadata', {}).get('pmid')}
                        evidencia_filtrada = []
                        for doc in evidence:
                            pmid = doc.get('metadata', {}).get('pmid')
                            if pmid is None or pmid not in pmids_existentes:
                                evidencia_filtrada.append(doc)
                                if pmid:
                                    pmids_existentes.add(pmid)
                        all_evidence.extend(evidencia_filtrada)
                        logger.info(f"Encontrada evidencia para consulta '{query}': {len(evidence)} documentos")
                except Exception as e:
                    logger.error(f"Error buscando evidencia para '{query}': {str(e)}")
                    continue
            
            # Si no hay evidencia semántica, intentar búsqueda general
            if not all_evidence:
                logger.info("No hay evidencia semántica, intentando búsqueda general")
                evidence = await self.semantic_search_service.search_medical_evidence(symptoms, max_results=5)
                if evidence:
                    all_evidence.extend(evidence)
            
            state["medical_evidence"] = all_evidence
            logger.info(f"Evidencia médica encontrada en iteración {state['search_iterations']}: {len(all_evidence)} documentos")
            
            return state
            
        except Exception as e:
            logger.error(f"Error en búsqueda semántica: {str(e)}")
            state["error"] = f"Error en búsqueda semántica: {str(e)}"
            return state

    async def _evaluate_relevance(self, state: DiagnosisState) -> DiagnosisState:
        """Evalúa la relevancia de la evidencia encontrada"""
        try:
            medical_evidence = state["medical_evidence"]
            symptoms = state["symptoms"]
            
            evaluation = await self.evaluation_tools.evaluate_relevance(symptoms, medical_evidence)
            state["evidence_evaluation"] = evaluation
            
            return state
            
        except Exception as e:
            logger.error(f"Error evaluando relevancia: {str(e)}")
            state["error"] = f"Error evaluando relevancia: {str(e)}"
            return state
    
    def _evaluate_relevance_router(self, state: DiagnosisState) -> str:
        """Determina el siguiente paso basado en la evaluación de relevancia"""
        try:
            evaluation = state.get("evidence_evaluation", {})
            iterations = state.get("search_iterations", 0)
            max_iterations = state.get("max_iterations", 3)
            
            return self.evaluation_tools.evaluate_relevance_router(evaluation, iterations, max_iterations)
            
        except Exception as e:
            logger.error(f"Error en lógica condicional: {str(e)}")
            return "expand"  # Por defecto, expandir
    
    async def _generate_pubmed_queries(self, state: DiagnosisState) -> DiagnosisState:
        """Genera consultas de búsqueda PubMed optimizadas para fisioterapia"""
        try:
            symptoms = state["symptoms"]
            age = state["patient_age"]
            gender = state["patient_gender"]
            
            queries = await self.search_tools.generate_pubmed_queries(symptoms, age, gender)
            state["pubmed_queries"] = queries
            
            if not queries:
                state["error"] = "Error generando consultas PubMed"
            
            return state

        except Exception as e:
            logger.error(f"Error generando consultas PubMed: {str(e)}")
            state["error"] = f"Error generando consultas PubMed: {str(e)}"
            return state
    
    async def _search_pubmed(self, state: DiagnosisState) -> DiagnosisState:
        """Realiza búsqueda en PubMed"""
        try:
            pubmed_queries = state["pubmed_queries"]
            all_evidence = []
            
    
            # Buscar evidencia en PubMed para cada consulta
            for query in pubmed_queries[:3]:  # Limitar a 3 consultas por iteración
                try:
                    evidence = await self.pubmed_service.search_and_fetch(query, max_results=3)
                    if evidence:
                        all_evidence.extend(evidence)
                        logger.info(f"Encontrada evidencia para consulta '{query}': {len(evidence)} documentos")
                except Exception as e:
                    logger.error(f"Error buscando evidencia para '{query}': {str(e)}")
                    continue
            
            logger.info(f"Evidencia médica encontrada en iteración {state['search_iterations']}: {len(all_evidence)} documentos")
            
            return state
            
        except Exception as e:
            logger.error(f"Error en búsqueda en PubMed: {str(e)}")
            state["error"] = f"Error en búsqueda en PubMed: {str(e)}"
            return state
    
    async def _generate_diagnosis(self, state: DiagnosisState) -> DiagnosisState:
        """Genera la evaluación completa del paciente para el fisioterapeuta"""
        try:
            symptoms = state["symptoms"]
            patient_age = state["patient_age"]
            patient_gender = state["patient_gender"]
            medical_evidence = state["medical_evidence"]
            evidence_evaluation = state.get("evidence_evaluation", {})
            
            patient_evaluation = await self.diagnosis_tools.generate_diagnosis(
                symptoms, patient_age, patient_gender, medical_evidence, evidence_evaluation
            )
            
            state["patient_evaluation"] = patient_evaluation
            return state
            
        except Exception as e:
            logger.error(f"Error generando evaluación del paciente: {str(e)}")
            state["error"] = f"Error generando evaluación del paciente: {str(e)}"
            return state
    
    async def _no_evidence_found(self, state: DiagnosisState) -> DiagnosisState:
        """Genera una respuesta cuando no se encuentra evidencia médica"""
        try:
            state["patient_evaluation"] = self.diagnosis_tools.generate_no_evidence_response()
            return state
        except Exception as e:
            logger.error(f"Error generando respuesta cuando no se encuentra evidencia médica: {str(e)}")
            state["error"] = f"Error generando respuesta cuando no se encuentra evidencia médica: {str(e)}"
            return state
    

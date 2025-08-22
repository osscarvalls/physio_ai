"""
Servicio para manejar operaciones con modelos de lenguaje
"""

import logging
from typing import List, Dict, Any
from langchain_openai import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.chains.combine_documents import create_stuff_documents_chain
from langchain.schema import Document

from app.config.settings import settings

logger = logging.getLogger(__name__)


class LLMService:
    """Servicio para manejar operaciones con LLMs"""
    
    def __init__(self):
        """Inicializa el servicio LLM"""
        self.llm = ChatOpenAI(
            model=settings.OPENAI_MODEL,
            temperature=settings.OPENAI_TEMPERATURE,
            api_key=settings.OPENAI_API_KEY
        )
        self._setup_chains()
    
    def _setup_chains(self):
        """Configura las cadenas de procesamiento"""
        self.document_chain = create_stuff_documents_chain(
            llm=self.llm,
            prompt=self._create_diagnosis_prompt()
        )
    
    def _create_diagnosis_prompt(self) -> ChatPromptTemplate:
        """Crea el prompt para diagnóstico médico"""
        return ChatPromptTemplate.from_messages([
            ("system", self._get_diagnosis_system_prompt()),
            ("human", "{input}")
        ])
    
    def _get_diagnosis_system_prompt(self) -> str:
        """Retorna el prompt del sistema para diagnóstico"""
        return """Eres un asistente médico experto. Tu tarea es:

1. Analizar los síntomas proporcionados por el paciente
2. Generar un diagnóstico preliminar basado en la evidencia médica disponible
3. Proporcionar recomendaciones médicas apropiadas
4. Siempre recordar que eres un asistente y no reemplazas la consulta médica profesional

IMPORTANTE: 
- Sé preciso pero conservador en tus diagnósticos
- Siempre recomienda consultar con un profesional médico
- Usa un lenguaje claro y comprensible
- Basa tus respuestas en la evidencia médica proporcionada

Contexto médico disponible:
{context}

Por favor, analiza los siguientes síntomas: {input}"""
    
    async def generate_diagnosis(
        self, 
        symptoms: str, 
        medical_context: List[Document]
    ) -> Dict[str, Any]:
        """
        Genera un diagnóstico basado en síntomas y contexto médico
        
        Args:
            symptoms: Descripción de los síntomas
            medical_context: Documentos médicos relevantes
            
        Returns:
            Dict con el diagnóstico y metadatos
        """
        try:
            if not symptoms or not symptoms.strip():
                return {
                    "diagnosis": "Por favor, proporcione síntomas para el diagnóstico",
                    "confidence": 0.0,
                    "error": "Síntomas no proporcionados"
                }
            
            # Procesar con la cadena de documentos
            result = await self.document_chain.ainvoke({
                "input": symptoms,
                "context": self._format_context(medical_context)
            })
            
            return {
                "diagnosis": result.get("answer", "No se pudo generar diagnóstico"),
                "confidence": 0.8,  # Placeholder - se puede implementar scoring real
                "timestamp": result.get("timestamp"),
                "context_used": len(medical_context)
            }
            
        except Exception as e:
            logger.error(f"Error generando diagnóstico: {str(e)}")
            return {
                "diagnosis": "No se pudo generar diagnóstico en este momento",
                "confidence": 0.0,
                "error": str(e)
            }
    
    def _format_context(self, documents: List[Document]) -> str:
        """Formatea el contexto médico para el prompt"""
        if not documents:
            return "No hay evidencia médica disponible."
        
        context_parts = []
        for i, doc in enumerate(documents[:5], 1):  # Limitar a 5 documentos
            context_parts.append(f"Documento {i}: {doc.page_content[:500]}...")
        
        return "\n\n".join(context_parts)
    
    async def generate_search_queries(self, symptoms: str) -> List[str]:
        """
        Genera consultas de búsqueda optimizadas para PubMed
        
        Args:
            symptoms: Síntomas del paciente
            
        Returns:
            Lista de consultas de búsqueda
        """
        try:
            search_prompt = ChatPromptTemplate.from_template("""
            Eres un experto en búsquedas médicas. Genera 3-5 consultas de búsqueda 
            optimizadas para PubMed basadas en los siguientes síntomas.
            
            Síntomas: {symptoms}
            
            Genera solo las consultas, una por línea, sin numeración ni formato adicional.
            """)
            
            chain = search_prompt | self.llm
            result = await chain.ainvoke({"symptoms": symptoms})
            
            # Parsear resultado
            queries = [
                query.strip() 
                for query in result.content.split('\n') 
                if query.strip()
            ]
            
            return queries[:5]  # Limitar a 5 consultas
            
        except Exception as e:
            logger.error(f"Error generando consultas de búsqueda: {str(e)}")
            return [symptoms]  # Fallback a síntomas originales

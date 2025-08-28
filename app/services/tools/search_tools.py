"""
Herramientas de búsqueda para el servicio de diagnóstico
"""

import logging
from typing import List, Dict, Any
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI

logger = logging.getLogger(__name__)


class SearchTools:
    """Herramientas para generar consultas y realizar búsquedas"""
    
    def __init__(self, llm: ChatOpenAI):
        self.llm = llm
    
    async def generate_semantic_queries(self, symptoms: str, age: int, gender: str) -> List[str]:
        """Genera consultas de búsqueda optimizadas para fisioterapia"""
        try:
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta clínico especializado en la búsqueda de literatura científica relevante para casos de fisioterapia.

            Paciente: {age} años, {gender}
            Síntomas principales: {symptoms}

            Tu tarea es generar entre 3 y 5 consultas de búsqueda semántica altamente optimizadas para Qdrant, enfocadas en encontrar artículos científicos, revisiones sistemáticas o guías clínicas que ayuden a:
            - Identificar diagnósticos diferenciales relevantes en fisioterapia para estos síntomas.
            - Proponer pruebas de evaluación clínica y funcional específicas.
            - Sugerir intervenciones y tratamientos fisioterapéuticos basados en evidencia.
            - Localizar evidencia clínica reciente y de alta calidad.

            Las consultas deben:
            - Ser claras, concisas y específicas, usando términos técnicos y palabras clave relevantes.
            - Incluir, si es posible, combinaciones de síntomas, contexto clínico y términos de fisioterapia.
            - Evitar frases genéricas o demasiado amplias.
            - No incluir numeración, formato adicional ni explicaciones. Escribe solo una consulta por línea.
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
            
            logger.info(f"Consultas de búsqueda semántica generadas: {len(queries)}")
            return queries[:5]  # Limitar a 5 consultas
            
        except Exception as e:
            logger.error(f"Error generando consultas semánticas: {str(e)}")
            return []
    
    async def generate_pubmed_queries(self, symptoms: str, age: int, gender: str) -> List[str]:
        """Genera consultas de búsqueda PubMed optimizadas para fisioterapia"""
        try:
            prompt = ChatPromptTemplate.from_template("""
            Eres un fisioterapeuta clínico especializado en la búsqueda de literatura científica relevante para casos de fisioterapia.

            Paciente: {age} años, {gender}
            Síntomas principales: {symptoms}
            
            Tu tarea es generar entre 3 y 5 consultas de búsqueda PubMed altamente optimizadas para encontrar artículos científicos, revisiones sistemáticas o guías clínicas que ayuden a:
            - Identificar diagnósticos diferenciales relevantes en fisioterapia para estos síntomas.
            - Proponer pruebas de evaluación clínica y funcional específicas.
            - Sugerir intervenciones y tratamientos fisioterapéuticos basados en evidencia.
            - Localizar evidencia clínica reciente y de alta calidad.
            
            Las consultas deben:
            - Ser claras, concisas y específicas, usando términos técnicos y palabras clave relevantes.
            - Incluir, si es posible, combinaciones de síntomas, contexto clínico y términos de fisioterapia.
            - Evitar frases genéricas o demasiado amplias.
            - No incluir numeración, formato adicional ni explicaciones. Escribe solo una consulta por línea.
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
            
            logger.info(f"Consultas de búsqueda PubMed generadas: {len(queries)}")
            return queries[:5]  # Limitar a 5 consultas

        except Exception as e:
            logger.error(f"Error generando consultas PubMed: {str(e)}")
            return []

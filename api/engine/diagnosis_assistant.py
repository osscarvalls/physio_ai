from api.engine.utils import get_prompts, set_openai
from api.engine.evidence_retrieval import EvidenceRetrievalEngine

from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
from langchain_openai import ChatOpenAI

class DiagnosisAssistantEngine():
    def __init__(self) -> None:
        set_openai()
        self.prompts = get_prompts()
        # Create the prompt template
        self.prompt = PromptTemplate(
            input_variables=["symptoms", "context"],
            template=self.prompts["diagnosis_assistant"]
        )

        # Create the RAG chain
        self.llm = ChatOpenAI(model='gpt-4o-mini',temperature=0)
        self.rag_chain = LLMChain(
            llm=self.llm,
            prompt=self.prompt,
            verbose=False
        )

        self.evidence_retrieval = EvidenceRetrievalEngine()

    def generate_diagnosis(self,symptoms):

        self.evidence_retrieval(symptoms)

        retrieved_docs = self.evidence_retrieval(symptoms)

        # Example usage of the RAG chain
        context = "\n".join([doc.page_content for doc in retrieved_docs])
        result = self.rag_chain.invoke({"symptoms": symptoms, "context": context})

        return result['text']
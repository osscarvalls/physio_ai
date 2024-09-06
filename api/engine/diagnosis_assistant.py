from api.engine.utils import get_prompts, set_openai
from api.engine.evidence_retrieval import EvidenceRetrievalEngine

from langchain.prompts import ChatPromptTemplate
from langchain.chains import create_retrieval_chain
from langchain_openai import ChatOpenAI
from langchain.chains.combine_documents import create_stuff_documents_chain

class DiagnosisAssistantEngine():
    def __init__(self) -> None:
        set_openai()
        self.system_prompt = get_prompts()
        # Create the prompt template
        self.prompt = ChatPromptTemplate.from_messages(
            [
                ("system", self.system_prompt["diagnosis_assistant"]),
                ("human", "{input}"),
            ]
        )

        # Define the language model
        self.llm = ChatOpenAI(model='gpt-4o-mini',temperature=0)

        # Create the question-answering chain
        self.chain = create_stuff_documents_chain(self.llm, self.prompt)

        # Create the evidence retrieval engine
        self.evidence_retrieval = EvidenceRetrievalEngine()

    def generate_diagnosis(self,symptoms):

        # Retrieve the evidence
        retriever = self.evidence_retrieval(symptoms)

        rag_chain = create_retrieval_chain(retriever, self.chain)

        # Example usage of the RAG chain
        result = rag_chain.invoke({"input": symptoms})
        print(result)

        return result['answer']
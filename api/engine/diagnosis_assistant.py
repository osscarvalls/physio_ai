from engine.utils import get_prompts, set_openai
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain


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
        self.rag_chain = LLMChain(
            llm=self.llm,
            prompt=self.prompt,
            verbose=False
        )

    def retrieve_evidence(self):
        retriever = vectorstore.as_retriever(search_type="similarity", search_kwargs={"k": 10})


    def generate_diagnosis(self,symptoms):

        retrieved_docs = self.retrieve_evidence(symptoms)

        # Example usage of the RAG chain
        context = "\n".join([doc.page_content for doc in retrieved_docs])
        result = self.rag_chain.invoke({"symptoms": symptoms, "context": context})

        return result['text']
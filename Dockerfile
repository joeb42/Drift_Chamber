FROM continuumio/miniconda3

WORKDIR usr/local/app

COPY . .
RUN conda install pip && pip install . 

CMD ["python", "main.py"]

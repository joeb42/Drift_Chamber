FROM continuumio/miniconda3

WORKDIR usr/local/app

COPY . .
#RUN conda init bash
RUN conda install pip && pip install -e . 

CMD ["python", "main.py"]

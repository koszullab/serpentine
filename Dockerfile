FROM python:3

WORKDIR /serpentine

ADD serpentine /serpentine
ADD requirements.txt /serpentine
ADD setup.py /serpentine

RUN pip install -Ur /serpentine

ENTRYPOINT ["serpentine"]
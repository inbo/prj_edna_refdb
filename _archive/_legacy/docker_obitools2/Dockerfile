FROM python:2

LABEL author="Roman Luštrik (roman.lustrik@biolitika.si)" #github.com/romunov/obitools
LABEL author="Pieter Verschelde (pieter.verschelde@inbo.be)"
# Install OBITools.
RUN wget https://pythonhosted.org/OBITools/_downloads/get-obitools.py && python get-obitools.py

RUN chmod -R 700 OBITools-1.2.13
RUN chmod 700 obitools

RUN mv obitools /usr/local/bin

RUN rm get-obitools.py

ENV PATH=/OBITools-1.2.13/bin:$PATH

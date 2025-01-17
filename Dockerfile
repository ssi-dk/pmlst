FROM python:3.12-slim-bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN apt update -qq; \
    apt install -y git \
    build-essential \
    apt-utils \
    wget \
    ncbi-blast+ \
    libz-dev \
    procps \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
    
RUN python -m pip install --upgrade pip

ENV DEBIAN_FRONTEND Teletype

# Install python dependencies
RUN pip install -U biopython==1.85 tabulate==0.9.0 cgecore==2.0.0;

# Install kma
RUN cd /usr/src; \
    git clone --depth 1 -b 1.4.15 https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma && make; \
    mv kma* /usr/bin/; \
    cd ..; \
    rm -rf kma/;

COPY pmlst.py /usr/src/pmlst.py 

RUN chmod 755 /usr/src/pmlst.py; 

# Install database
RUN cd /; \
    mkdir databases; \
    cd /databases/; \
    git clone https://bitbucket.org/genomicepidemiology/pmlst_db.git; \
    cd /databases/pmlst_db; \
    python INSTALL.py; \
    rm -rf .git;

# Environmental variables
ENV PMLST_DB /databases/pmlst_db/
RUN echo "Database path is $PMLST_DB"

ENV PATH $PATH:/usr/src
# Setup .bashrc file for convenience during debugging
RUN echo "alias ls='ls -h --color=tty'\n"\
"alias ll='ls -lrt'\n"\
"alias l='less'\n"\
"alias du='du -hP --max-depth=1'\n"\
"alias cwd='readlink -f .'\n"\
"PATH=$PATH\n">> ~/.bashrc

# Setup environment
RUN cd /; \
    mkdir app;
WORKDIR /app

# Execute program when running the container
ENTRYPOINT ["/usr/src/pmlst.py"]

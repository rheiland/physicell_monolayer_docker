# Use Ubuntu as base image
FROM ubuntu:22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install build essentials
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app
#RUN mkdir -p /opt/pc_monolayer/config &&\
#    chmod -R 777 /opt/pc_monolayer

# Copy source files
COPY . .
#COPY ./config/* /opt/pc_monolayer/config/
#COPY ./rwh2.xml .
#COPY ./config/* .
#COPY ./rwh2.xml ./config/PhysiCell.xml

# Compile the project
RUN make
#COPY ./project /opt/pc_monolayer/
#COPY ./project .
RUN chmod +x project


# Default command (run the compiled program or bash)
#CMD ["make"]
#CMD ["./project ./rwh2.xml"]
#CMD ["./project rwh2.xml"]
CMD ["./project"]
#ENTRYPOINT ["./project"]

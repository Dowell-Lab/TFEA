DOCKER_BIND_TEST_PATH=$(CURDIR)/docker_test_out
SINGULARITY_BIND_TEST_PATH=$(CURDIR)/singularity_test_out
IMAGE_TEST_PATH=/usr/local/lib/python3.8/dist-packages/TFEA/test/test_files/test_output

.PHONY: all test docker docker_test singularity singularity_test

all: singularity docker

test: docker_test singularity_test

docker:
	docker build -t jdrubin/tfea -f $(CURDIR)/Docker/Dockerfile $(CURDIR)

docker_test: docker
	mkdir -p $(DOCKER_BIND_TEST_PATH)
	docker run --rm --mount type=bind,source=$(DOCKER_BIND_TEST_PATH),target=$(IMAGE_TEST_PATH) jdrubin/tfea TFEA --test-full

singularity: tfea.sif

tfea.sif:
	singularity build -f tfea.sif Singularity/TFEA.def

singularity_test: singularity
	mkdir -p $(SINGULARITY_BIND_TEST_PATH)
	singularity exec --bind $(SINGULARITY_BIND_TEST_PATH):$(IMAGE_TEST_PATH) tfea.sif TFEA --test-full

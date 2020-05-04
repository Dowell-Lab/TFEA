BIND_TEST_PATH=$(CURDIR)/test_out
IMAGE_TEST_PATH=/usr/local/lib/python3.6/dist-packages/TFEA/test/test_files/test_output/

.PHONY: all test docker docker_test singularity singularity_test

all: singularity docker

test: docker_test singularity_test

docker:
	docker build -t jdrubin/tfea -f $(CURDIR)/Docker/Dockerfile $(CURDIR)

docker_test: docker
	docker run --rm --mount type=bind,source=$(BIND_TEST_PATH),target=$(IMAGE_TEST_PATH) jdrubin/tfea TFEA --test-full

singularity:
	singularity build tfea.sif Singularity/TFEA.def

singularity_test: singularity
	singularity exec --bind $(BIND_TEST_PATH):$(IMAGE_TEST_PATH) tfea.sif TFEA --test-full

COURSE=DNPM_2023_autumn
CODE?=2_all_velocities
SCRIPT?=task2

compile:
	g++ src/$(COURSE)/$(CODE).cpp -o bin/$(CODE)

run: compile
	bin/$(CODE)

gif: run
	gnuplot -c  gnuplot/$(SCRIPT).gnu

mp4: gif
	ffmpeg -y -i data/$(SCRIPT).gif -movflags faststart -pix_fmt yuv420p data/$(SCRIPT).mp4

all: mp4

clean:
	echo "Cleaned!"


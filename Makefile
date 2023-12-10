COURSE=DNPM_2023_autumn
CODE?=3_2d_peak
SCRIPT?=task3

compile:
	g++ src/$(COURSE)/$(CODE).cpp -o bin/$(CODE)

run: compile
	bin/$(CODE)

gif: run
	gnuplot -c  gnuplot/$(SCRIPT).gnu

mp4: gif
	ffmpeg -y -i data/$(SCRIPT).gif -movflags faststart -pix_fmt yuv420p data/$(SCRIPT).mp4

total: 
	gnuplot -c gnuplot/plain.gnu

all: mp4 total

clean:
	echo "Cleaned!"


package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"math/rand"
)

func SimulateEvolutionFourByFour(mutationMatrix [][]float64, seqA fasta.Fasta) fasta.Fasta {
	var answer fasta.Fasta = fasta.Fasta{Name: "seqB", Seq: make([]dna.Base, len(seqA.Seq))}
	dna.AllToUpper(seqA.Seq)
	for currPos := range seqA.Seq {

		substitutedBase := substitute(seqA.Seq[currPos], mutationMatrix)
		answer.Seq[currPos] = substitutedBase
	}
	return answer

}
func substitute(inBase dna.Base, mutationMatrix [][]float64) dna.Base {
	var currRand float64 = rand.Float64()
	if inBase > 3 {

		return inBase
	}
	if currRand < mutationMatrix[inBase][dna.A] {
		return dna.A
	} else if currRand < mutationMatrix[inBase][dna.A]+mutationMatrix[inBase][dna.C] {
		return dna.C
	} else if currRand < mutationMatrix[inBase][dna.A]+mutationMatrix[inBase][dna.C]+mutationMatrix[inBase][dna.G] {
		return dna.G
	}
	return dna.T
}

func main() {
	var mySeq fasta.Fasta = fasta.Fasta{Name: "mySeq", Seq: dna.StringToBases("GG--CATCCGTnnnnCACaaACCggtAT-G")}
	var mutMatrix [][]float64 = [][]float64{[]float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}}
	fmt.Println(dna.BasesToString(mySeq.Seq))
	answer := SimulateEvolutionFourByFour(mutMatrix, mySeq)
	fmt.Println(dna.BasesToString(answer.Seq))

}

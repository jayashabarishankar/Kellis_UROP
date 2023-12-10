package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var SimulateEvolutionFourByFourTest = []struct {
	InSeq          fasta.Fasta
	ExpectedSeq    fasta.Fasta
	MutationMatrix [][]float64
}{
	{InSeq: fasta.Fasta{Name: "mySeq", Seq: dna.StringToBases("GG--CATCCGTnnnnCACaaACCggtAT-G")},
		ExpectedSeq:    fasta.Fasta{Name: "seqB", Seq: dna.StringToBases("CC--CCCCCCCNNNNCCCCCCCCCCCCC-C")},
		MutationMatrix: [][]float64{[]float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}, []float64{0, 1, 0, 0}},
	},
}

func TestSimulateEvolutionFourByFour(t *testing.T) {
	var ObservedSeq fasta.Fasta
	for _, v := range SimulateEvolutionFourByFourTest {
		ObservedSeq = SimulateEvolutionFourByFour(v.MutationMatrix, v.InSeq)
		if !fasta.IsEqual(ObservedSeq, v.ExpectedSeq) {
			t.Errorf("Error: output of SimulateEvolutionFourByFour was not expected.\n")
		}

	}
}

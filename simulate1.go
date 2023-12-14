package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"gonum.org/v1/gonum/mat"
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

type SubstitutionTree struct {
	NodeName       string
	BranchLength   float64     //this is the evolutionary time distance between this node and Up
	MutationMatrix [][]float64 //this is the substitution probability matrix between this node and Up
	Left           *SubstitutionTree
	Right          *SubstitutionTree
	Up             *SubstitutionTree
}

func BuildSubstitutionTree(startingMatrix [][]float64, startingBranchLength float64, root *SubstitutionTree) {
	//var matrixExp float64 = root.Left.BranchLength / startingBranchLength

}

func floatMatrixToDenseMatrix(inMat [][]float64) *mat.Dense {
	rows := len(inMat)
	if rows < 1 {
		return mat.NewDense(0, 0, nil)
	}
	cols := len(inMat[0])
	dense := mat.NewDense(rows, cols, nil)
	for currRow := 0; currRow < rows; currRow++ {
		for currCol := 0; currCol < cols; currCol++ {
			dense.Set(currRow, currCol, inMat[currRow][currCol])
		}
	}
	return dense
}
func denseMatrixToFloatMatrix(dense *mat.Dense) [][]float64 {
	rows, cols := dense.Dims()
	var currRow, currCol int

	// Create a new float matrix with the same dimensions as the dense matrix
	floatMat := make([][]float64, rows)
	for i := range floatMat {
		floatMat[i] = make([]float64, cols)
	}

	// Populate the float matrix with values from the dense matrix
	for currRow = 0; currRow < rows; currRow++ {
		for currCol = 0; currCol < cols; currCol++ {
			floatMat[currRow][currCol] = dense.At(currRow, currCol)
		}
	}

	return floatMat
}
func main() {
	//var mySeq fasta.Fasta = fasta.Fasta{Name: "mySeq", Seq: dna.StringToBases("GG--CATCCGTnnnnCACaaACCggtAT-G")}
	var mutMatrix [][]float64 = [][]float64{[]float64{0.97, 0.01, 0.01, 0.01}, []float64{0.01, 0.97, 0.01, 0.01}, []float64{0.01, 0.01, 0.97, 0.01}, []float64{0.01, 0.01, 0.01, .97}}
	//fmt.Println(dna.BasesToString(mySeq.Seq))
	//answer := SimulateEvolutionFourByFour(mutMatrix, mySeq)
	//fmt.Println(dna.BasesToString(answer.Seq))
	myMat := floatMatrixToDenseMatrix(mutMatrix)
	myMat.Scale(2, myMat)
	myMat.Exp(myMat)
	//fmt.Println(mat.Formatted(myMat))
	myVec := floatMatrixToDenseMatrix([][]float64{{1}, {0}, {0}, {0}})
	var myAnswer mat.Dense
	myAnswer.Mul(myMat, myVec)
	fmt.Println(denseMatrixToFloatMatrix(&myAnswer))
	//myMat1 := denseMatrixToFloatMatrix(myMat)
	//fmt.Println(myMat1)
}

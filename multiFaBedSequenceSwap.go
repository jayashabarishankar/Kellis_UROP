package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
)

// MultiFa represents a collection of sequences in a single file (similar to a Fasta, but it is from different species)
type MultiFa struct {
	FilePath string
}

// SubsequenceSwap swaps a specified region between two sequences in a MultiFa object.
func SubsequenceSwap(mfa *MultiFa, background, foreground int, swapRegion bed.Bed) (*fasta.Fasta, error) {
	// Load the original sequences from the multiFa file.
	originalSeqs := fasta.Read(mfa.FilePath) //Converting a multiFa to a Fasta file

	// Ensure the specified sequences are within the valid range.
	if background < 0 || background >= len(originalSeqs) || foreground < 0 || foreground >= len(originalSeqs) {
		return nil, fmt.Errorf("Invalid sequence index")
	}

	// Swap the specified region between background and foreground sequences.
	backgroundSeq := originalSeqs[background].Seq
	foregroundSeq := originalSeqs[foreground].Seq

	start := swapRegion.ChromStart
	end := swapRegion.ChromEnd
	if start < 0 || end > len(backgroundSeq) || start >= end {
		return nil, fmt.Errorf("Invalid swap region")
	}

	// Perform the subsequence swap.
	backgroundSeq = append(backgroundSeq[:start], append(foregroundSeq[start:end], backgroundSeq[end:]...)...)
	foregroundSeq = append(foregroundSeq[:start], append(backgroundSeq[start:end], foregroundSeq[end:]...)...)

	originalSeqs[background].Seq = backgroundSeq
	originalSeqs[foreground].Seq = foregroundSeq

	// Create and return the new Fasta.
	return &fasta.Fasta{
		Seq: originalSeqs,
	}, nil
}

func main() {
	// Example usage of SubsequenceSwap function
	// Replace the arguments with actual values

	// Example multiFa
	mfa := &MultiFa{
		FilePath: "gonomics/cmd/multFaVisualizer/testdata/test.fa",
	}

	// Example swap region
	swapRegion := bed.Bed{
		ChromStart: 1,
		ChromEnd:   3,
	}

	// Specify the background and foreground sequences (adjust these indices accordingly)
	backgroundSeqIndex := 0
	foregroundSeqIndex := 1

	result, err := SubsequenceSwap(mfa, backgroundSeqIndex, foregroundSeqIndex, swapRegion)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}

	// Do something with the result, e.g., print it or further processing
	fmt.Println(result)
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

// Constants
#define MAX_SEQ_LENGTH 1000
#define MAX_SEQS 100
#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2
#define MIN_OVERLAP 20
#define OVERLAP_THRESHOLD 0.8

// Data structures
typedef struct {
    char id[50];
    char sequence[MAX_SEQ_LENGTH];
    int length;
} DNASequence;

typedef struct {
    int score;
    char alignmentA[MAX_SEQ_LENGTH * 2];
    char alignmentB[MAX_SEQ_LENGTH * 2];
    int alignmentLength;
} AlignmentResult;

typedef struct {
    int seq1_idx;
    int seq2_idx;
    int overlap_score;
    int overlap_length;
    int direction; // 0 = no overlap, 1 = seq1 before seq2, 2 = seq2 before seq1
} Overlap;

// Global variables
DNASequence sequences[MAX_SEQS];
int num_sequences = 0;
Overlap overlaps[MAX_SEQS][MAX_SEQS];

// Utility functions
void readSequencesFromFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(1);
    }

    char line[MAX_SEQ_LENGTH];
    while (fgets(line, sizeof(line), file) && num_sequences < MAX_SEQS) {
        
        // Remove newline character
        line[strcspn(line, "\n")] = '\0';
        
        // Check if line starts with '>' (header)
        if (line[0] == '>') {
            strncpy(sequences[num_sequences].id, line + 1, sizeof(sequences[num_sequences].id) - 1);
        } else {
            strncpy(sequences[num_sequences].sequence, line, sizeof(sequences[num_sequences].sequence) - 1);
            sequences[num_sequences].length = strlen(sequences[num_sequences].sequence);
            num_sequences++;
        }
    }

    fclose(file);
}

void writeResultsToFile(const char* filename, const char* content) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }
    fprintf(file, "%s", content);
    fclose(file);
}

void printSequence(DNASequence seq) {
    printf("ID: %s\n", seq.id);
    printf("Length: %d\n", seq.length);
    printf("Sequence: %.50s...\n", seq.sequence);
}

// Scoring function
int score(char a, char b) {
    return (a == b) ? MATCH_SCORE : MISMATCH_SCORE;
}

// Sequence alignment algorithms
AlignmentResult needlemanWunsch(const char* seqA, const char* seqB) {
    int lenA = strlen(seqA);
    int lenB = strlen(seqB);
    
    // Allocate DP table
    int** dp = (int**)malloc((lenA + 1) * sizeof(int*));
    for (int i = 0; i <= lenA; i++) {
        dp[i] = (int*)malloc((lenB + 1) * sizeof(int));
    }
    
    // Initialize DP table
    for (int i = 0; i <= lenA; i++) {
        dp[i][0] = i * GAP_PENALTY;
    }
    for (int j = 0; j <= lenB; j++) {
        dp[0][j] = j * GAP_PENALTY;
    }
    
    // Fill DP table
    for (int i = 1; i <= lenA; i++) {
        for (int j = 1; j <= lenB; j++) {
            int match = dp[i-1][j-1] + score(seqA[i-1], seqB[j-1]);
            int delete = dp[i-1][j] + GAP_PENALTY;
            int insert = dp[i][j-1] + GAP_PENALTY;
            dp[i][j] = (match > delete) ? match : delete;
            dp[i][j] = (dp[i][j] > insert) ? dp[i][j] : insert;
        }
    }
    
    // Backtrack to find alignment
    AlignmentResult result;
    result.score = dp[lenA][lenB];
    result.alignmentLength = 0;
    
    int i = lenA, j = lenB;
    int a_pos = 0, b_pos = 0;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + score(seqA[i-1], seqB[j-1])) {
            result.alignmentA[a_pos++] = seqA[i-1];
            result.alignmentB[b_pos++] = seqB[j-1];
            i--;
            j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + GAP_PENALTY) {
            result.alignmentA[a_pos++] = seqA[i-1];
            result.alignmentB[b_pos++] = '-';
            i--;
        } else {
            result.alignmentA[a_pos++] = '-';
            result.alignmentB[b_pos++] = seqB[j-1];
            j--;
        }
        result.alignmentLength++;
    }
    
    // Reverse the alignments
    for (int k = 0; k < result.alignmentLength / 2; k++) {
        char temp = result.alignmentA[k];
        result.alignmentA[k] = result.alignmentA[result.alignmentLength - 1 - k];
        result.alignmentA[result.alignmentLength - 1 - k] = temp;
        
        temp = result.alignmentB[k];
        result.alignmentB[k] = result.alignmentB[result.alignmentLength - 1 - k];
        result.alignmentB[result.alignmentLength - 1 - k] = temp;
    }
    
    result.alignmentA[result.alignmentLength] = '\0';
    result.alignmentB[result.alignmentLength] = '\0';
    
    // Free DP table
    for (int i = 0; i <= lenA; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return result;
}

AlignmentResult smithWaterman(const char* seqA, const char* seqB) {
    int lenA = strlen(seqA);
    int lenB = strlen(seqB);
    
    // Allocate DP table
    int** dp = (int**)malloc((lenA + 1) * sizeof(int*));
    for (int i = 0; i <= lenA; i++) {
        dp[i] = (int*)malloc((lenB + 1) * sizeof(int));
    }
    
    // Initialize DP table
    for (int i = 0; i <= lenA; i++) {
        dp[i][0] = 0;
    }
    for (int j = 0; j <= lenB; j++) {
        dp[0][j] = 0;
    }
    
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    // Fill DP table
    for (int i = 1; i <= lenA; i++) {
        for (int j = 1; j <= lenB; j++) {
            int match = dp[i-1][j-1] + score(seqA[i-1], seqB[j-1]);
            int delete = dp[i-1][j] + GAP_PENALTY;
            int insert = dp[i][j-1] + GAP_PENALTY;
            dp[i][j] = (match > delete) ? match : delete;
            dp[i][j] = (dp[i][j] > insert) ? dp[i][j] : insert;
            dp[i][j] = (dp[i][j] > 0) ? dp[i][j] : 0;
            
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Backtrack to find alignment
    AlignmentResult result;
    result.score = max_score;
    result.alignmentLength = 0;
    
    int i = max_i, j = max_j;
    int a_pos = 0, b_pos = 0;
    
    while (dp[i][j] != 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + score(seqA[i-1], seqB[j-1])) {
            result.alignmentA[a_pos++] = seqA[i-1];
            result.alignmentB[b_pos++] = seqB[j-1];
            i--;
            j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + GAP_PENALTY) {
            result.alignmentA[a_pos++] = seqA[i-1];
            result.alignmentB[b_pos++] = '-';
            i--;
        } else {
            result.alignmentA[a_pos++] = '-';
            result.alignmentB[b_pos++] = seqB[j-1];
            j--;
        }
        result.alignmentLength++;
    }
    
    // Reverse the alignments
    for (int k = 0; k < result.alignmentLength / 2; k++) {
        char temp = result.alignmentA[k];
        result.alignmentA[k] = result.alignmentA[result.alignmentLength - 1 - k];
        result.alignmentA[result.alignmentLength - 1 - k] = temp;
        
        temp = result.alignmentB[k];
        result.alignmentB[k] = result.alignmentB[result.alignmentLength - 1 - k];
        result.alignmentB[result.alignmentLength - 1 - k] = temp;
    }
    
    result.alignmentA[result.alignmentLength] = '\0';
    result.alignmentB[result.alignmentLength] = '\0';
    
    // Free DP table
    for (int i = 0; i <= lenA; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return result;
}

// Overlap detection for genome assembly
void computeOverlaps() {
    for (int i = 0; i < num_sequences; i++) {
        for (int j = 0; j < num_sequences; j++) {
            if (i == j) {
                overlaps[i][j].overlap_score = 0;
                overlaps[i][j].direction = 0;
                continue;
            }
            
            // Try both directions
            AlignmentResult result1 = smithWaterman(sequences[i].sequence, sequences[j].sequence);
            AlignmentResult result2 = smithWaterman(sequences[j].sequence, sequences[i].sequence);
            
            if (result1.score > result2.score) {
                overlaps[i][j].overlap_score = result1.score;
                overlaps[i][j].overlap_length = result1.alignmentLength;
                overlaps[i][j].direction = 1; // i before j
            } else {
                overlaps[i][j].overlap_score = result2.score;
                overlaps[i][j].overlap_length = result2.alignmentLength;
                overlaps[i][j].direction = 2; // j before i
            }
            
            overlaps[i][j].seq1_idx = i;
            overlaps[i][j].seq2_idx = j;
        }
    }
}

// Genome assembly using OLC approach
char* assembleGenome() {
    // Simple greedy assembly for demonstration
    // In a real implementation, we'd use a more sophisticated approach
    
    bool used[MAX_SEQS] = {false};
    char* assembled = (char*)malloc(MAX_SEQ_LENGTH * MAX_SEQS * sizeof(char));
    assembled[0] = '\0';
    
    // Start with the longest sequence
    int current_seq = 0;
    for (int i = 1; i < num_sequences; i++) {
        if (sequences[i].length > sequences[current_seq].length) {
            current_seq = i;
        }
    }
    strcpy(assembled, sequences[current_seq].sequence);
    used[current_seq] = true;
    
    while (true) {
        int best_overlap = -1;
        int best_seq = -1;
        int best_direction = 0;
        int best_pos = 0;
        
        // Find the best unused sequence with maximum overlap
        for (int i = 0; i < num_sequences; i++) {
            if (!used[i]) {
                if (overlaps[current_seq][i].overlap_score > best_overlap) {
                    best_overlap = overlaps[current_seq][i].overlap_score;
                    best_seq = i;
                    best_direction = overlaps[current_seq][i].direction;
                }
            }
        }
        
        if (best_overlap < MIN_OVERLAP) break;
        
        // Merge the sequences based on overlap direction
        if (best_direction == 1) { // current_seq before best_seq
            int overlap_len = overlaps[current_seq][best_seq].overlap_length;
            strcat(assembled, sequences[best_seq].sequence + overlap_len);
        } else { // best_seq before current_seq
            int overlap_len = overlaps[best_seq][current_seq].overlap_length;
            char temp[MAX_SEQ_LENGTH * MAX_SEQS];
            strcpy(temp, sequences[best_seq].sequence);
            strcat(temp, assembled + overlap_len);
            strcpy(assembled, temp);
        }
        
        used[best_seq] = true;
        current_seq = best_seq;
    }
    
    return assembled;
}

// Mutation detection
void detectMutations(const char* reference, const char* sample, char* result) {
    AlignmentResult alignment = needlemanWunsch(reference, sample);
    
    int pos = 0;
    result[0] = '\0';
    
    for (int i = 0; i < alignment.alignmentLength; i++) {
        if (alignment.alignmentA[i] != alignment.alignmentB[i]) {
            if (alignment.alignmentA[i] == '-') {
                // Insertion
                pos += sprintf(result + pos, "Insertion at position %d: inserted '%c'\n", 
                              i, alignment.alignmentB[i]);
            } else if (alignment.alignmentB[i] == '-') {
                // Deletion
                pos += sprintf(result + pos, "Deletion at position %d: deleted '%c'\n", 
                              i, alignment.alignmentA[i]);
            } else {
                // Substitution
                pos += sprintf(result + pos, "Substitution at position %d: '%c' -> '%c'\n", 
                              i, alignment.alignmentA[i], alignment.alignmentB[i]);
            }
        }
    }
}

// Main function
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }
    
    // Read input sequences
    readSequencesFromFile(argv[1]);
    
    printf("Read %d sequences:\n", num_sequences);
    for (int i = 0; i < num_sequences; i++) {
        printf("Sequence %d: ", i+1);
        printSequence(sequences[i]);
    }
    
    // Compute overlaps between sequences
    printf("\nComputing overlaps...\n");
    clock_t start = clock();
    computeOverlaps();
    clock_t end = clock();
    printf("Overlap computation took %.2f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    
    // Assemble genome
    printf("\nAssembling genome...\n");
    char* assembled_genome = assembleGenome();
    printf("Assembled genome length: %zu\n", strlen(assembled_genome));
    printf("First 100 bases: %.100s...\n", assembled_genome);
    
    // Write assembled genome to file
    writeResultsToFile("assembled_genome.txt", assembled_genome);
    
    // Detect mutations (using first two sequences as example)
    if (num_sequences >= 2) {
        printf("\nDetecting mutations between sequence 1 and 2...\n");
        char mutation_results[5000];
        detectMutations(sequences[0].sequence, sequences[1].sequence, mutation_results);
        printf("%s", mutation_results);
        writeResultsToFile("mutations.txt", mutation_results);
    }
    
    // Free memory
    free(assembled_genome);
    
    return 0;
}



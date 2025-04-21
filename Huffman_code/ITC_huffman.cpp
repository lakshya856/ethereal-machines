#include <iostream>
#include <queue>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

struct Node {
    char symbol;
    double probability;
    vector<Node*> children;
    string code;
};

struct Compare {
    bool operator()(Node* a, Node* b) {
        return a->probability > b->probability;
    }
};

void generateCodes(Node* root, string code) {
    if (!root) return;
    if (root->symbol != '\0') {
        root->code = code;
        return;
    }
    for (size_t i = 0; i < root->children.size(); i++) {
        generateCodes(root->children[i], code + to_string(i));
    }
}

Node* buildHuffmanTree(vector<pair<char, double>> &symbols, int base) {
    priority_queue<Node*, vector<Node*>, Compare> pq;
    for (auto &s : symbols) {
        Node* newNode = new Node{s.first, s.second, {}};
        pq.push(newNode);
    }
    
    while (pq.size() > 1) {
        vector<Node*> group;
        for (int i = 0; i < base && !pq.empty(); i++) {
            group.push_back(pq.top());
            pq.pop();
        }
        
        Node* parent = new Node{'\0', 0, group};
        for (auto child : group) {
            parent->probability += child->probability;
        }
        pq.push(parent);
    }
    return pq.top();
}

void collectCodes(Node* root, map<char, string>& codes) {
    if (!root) return;
    if (root->symbol != '\0') {
        codes[root->symbol] = root->code;
    }
    for (auto child : root->children) {
        collectCodes(child, codes);
    }
}

double calculateEntropy(vector<pair<char, double>>& symbols, int base) {
    double entropy = 0;
    for (auto &s : symbols) {
        entropy -= s.second * (log(s.second) / log(base));
    }
    return entropy;
}

double averageCodeLength(map<char, string>& codes, vector<pair<char, double>>& symbols) {
    double avgLength = 0;
    for (auto &s : symbols) {
        avgLength += s.second * codes[s.first].length();
    }
    return avgLength;
}

string decode(Node* root, string encoded, map<char, string>& codes) {
    string decoded = "";
    Node* current = root;
    for (char c : encoded) {
        if (!current->children.empty()) {
            current = current->children[c - '0'];
        }
        if (current->symbol != '\0') {
            decoded += current->symbol;
            current = root;
        }
    }
    
    cout << "Decoded Huffman Code: ";
    for (char c : decoded) {
        cout << codes[c] << " ";
    }
    cout << endl;
    
    return decoded;
}

string encode(string message, map<char, string>& codes) {
    string encoded = "";
    for (char c : message) {
        encoded += codes[c];
    }
    return encoded;
}

int main() {
    int numSymbols;
    cout << "Enter number of symbols (max 10): ";
    cin >> numSymbols;
    if (numSymbols < 1 || numSymbols > 10) {
        cout << "Invalid number! Defaulting to 7 symbols." << endl;
        numSymbols = 7;
    }
    
    vector<pair<char, double>> symbols;
    double cumulativeProbability = 0.0;
    
    for (int i = 0; i < numSymbols; i++) {
        char symbol = 'A' + i;
        double probability;
        cout << "Enter probability for " << symbol << ": ";
        cin >> probability;
        cumulativeProbability += probability;
        if (cumulativeProbability > 1.0) {
            cout << "Total probability exceeds 1. Please re-enter probability for " << symbol << "!" << endl;
            cumulativeProbability -= probability;
            i--;
            continue;
        }
        symbols.push_back({symbol, probability});
        cout << "Probability left to be assigned for remaining variables is: " << 1 - cumulativeProbability << endl;
    }
    
    if (cumulativeProbability != 1.0) {
        cout << "Error: Sum of all probabilities must be exactly 1. Please restart and enter valid probabilities." << endl;
        return 1;
    }
    
    sort(symbols.begin(), symbols.end(), [](const pair<char, double>& a, const pair<char, double>& b) {
        return a.second > b.second;
    });
    
    int base;
    cout << "Enter base for Huffman coding: ";
    cin >> base;
    
    if (base < 2) {
        cout << "Invalid base! Defaulting to base-2." << endl;
        base = 2;
    }
    
    Node* root = buildHuffmanTree(symbols, base);
    generateCodes(root, "");
    
    map<char, string> codes;
    collectCodes(root, codes);
    
    cout << "\nGenerated Huffman Codes:" << endl;
    for (auto &entry : codes) {
        cout << entry.first << " : " << entry.second << endl;
    }
    
    double entropy = calculateEntropy(symbols, base);
    double avgLength = averageCodeLength(codes, symbols);
    double efficiency = (entropy / avgLength) * 100;
    
    cout << "\nEntropy: " << entropy << " bits/symbol" << endl;
    cout << "Average Code Length: " << avgLength << " bits" << endl;
    cout << "Coding Efficiency: " << (efficiency > 100 ? 100 : efficiency) << " %" << endl;
    
    string message;
    cout << "\nEnter a message to encode: ";
    cin >> message;
    string encodedMessage = encode(message, codes);
    cout << "Encoded Sequence: " << encodedMessage << endl;
    
    string encoded;
    cout << "\nEnter encoded sequence to decode: ";
    cin >> encoded;
    
    string decoded = decode(root, encoded, codes);
    cout << "Decoded Data: " << decoded << endl;
    
    return 0;
}

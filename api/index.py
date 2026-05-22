from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import requests
import re
import math
import time
from unidecode import unidecode
from datetime import datetime, timedelta
from urllib.parse import quote
import sys
import os
from collections import defaultdict

# Adiciona o diretório pai ao path para importar os módulos
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Importa os dicionários dos arquivos separados
from dicionario_quimico import DICIONARIO_QUIMICO
from compostos_especiais import COMPOSTOS_ESPECIAIS, FORMULAS_GRANDES

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": "*"}})

# Cache em memória
class Cache:
    def __init__(self, timeout_seconds=3600):
        self.data = {}
        self.timeout = timedelta(seconds=timeout_seconds)
    
    def get(self, key):
        if key in self.data:
            item, timestamp = self.data[key]
            if datetime.now() - timestamp < self.timeout:
                return item
            del self.data[key]
        return None
    
    def set(self, key, value):
        self.data[key] = (value, datetime.now())

cache_store = Cache(timeout_seconds=3600)

# ============================================
# TRADUTOR QUÍMICO
# ============================================
class TradutorQuimico:
    @staticmethod
    def traduzir(nome):
        nome_lower = nome.lower().strip()
        nome_sem_acento = unidecode(nome_lower)
        
        if nome_lower in DICIONARIO_QUIMICO:
            return DICIONARIO_QUIMICO[nome_lower]
        
        if nome_sem_acento in DICIONARIO_QUIMICO:
            return DICIONARIO_QUIMICO[nome_sem_acento]
        
        palavras = nome_lower.split()
        palavras_traduzidas = []
        for palavra in palavras:
            if palavra in DICIONARIO_QUIMICO:
                palavras_traduzidas.append(DICIONARIO_QUIMICO[palavra])
            elif palavra not in ['de', 'da', 'do', 'das', 'dos', 'e', 'a', 'o', 'as', 'os']:
                palavras_traduzidas.append(palavra)
        
        if palavras_traduzidas:
            traducao = ' '.join(palavras_traduzidas)
            if traducao != nome_lower:
                return traducao
        
        try:
            url = "https://translate.googleapis.com/translate_a/single"
            params = {
                'client': 'gtx',
                'sl': 'pt',
                'tl': 'en',
                'dt': 't',
                'q': nome
            }
            response = requests.get(url, params=params, timeout=2)
            if response.status_code == 200:
                data = response.json()
                return data[0][0][0]
        except:
            pass
        
        return nome

# ============================================
# DETECTOR DE TIPO
# ============================================
class DetectorTipoBusca:
    @staticmethod
    def detectar_tipo(query):
        query = query.strip()
        
        if query.isdigit():
            return 'CID'
        elif query.upper() in FORMULAS_GRANDES:
            return 'FORMULA_GRANDE'
        elif re.match(r'^[A-Za-z][A-Za-z0-9]*$', query):
            if len(re.findall(r'\d+', query)) > 5:
                return 'FORMULA_GRANDE'
            return 'FORMULA_SIMPLES'
        elif re.search(r'[\(\)\[\]·]', query):
            return 'FORMULA_COMPLEXA'
        elif ' ' in query or len(query) > 3:
            return 'NOME_QUIMICO'
        else:
            return 'DESCONHECIDO'

# ============================================
# API PUBCHEM
# ============================================
class PubChemAPI:
    
    @staticmethod
    def validar_ligacoes(atoms, bonds):
        """Valida e corrige ligações químicas (evita H com múltiplas ligações)"""
        if not bonds or len(bonds) == 0:
            return bonds
        
        # Contar ligações por átomo
        bond_count = defaultdict(int)
        for bond in bonds:
            bond_count[bond['atom1']] += 1
            bond_count[bond['atom2']] += 1
        
        # Valências máximas por elemento
        valencias = {
            'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'Cl': 1,
            'Br': 1, 'I': 1, 'S': 6, 'P': 5, 'B': 3, 'Si': 4,
            'Na': 1, 'Mg': 2, 'Ca': 2, 'K': 1, 'Fe': 6, 'Cu': 4,
            'Zn': 2, 'Ag': 1, 'Au': 3, 'Hg': 2, 'Pb': 4
        }
        
        # Identificar ligações problemáticas
        bonds_to_remove = set()
        
        for bond in bonds:
            atom1 = atoms[bond['atom1']]
            atom2 = atoms[bond['atom2']]
            
            # Regra 1: Hidrogênio só pode ter 1 ligação
            if atom1['element'] == 'H' and bond_count[bond['atom1']] > 1:
                bonds_to_remove.add((min(bond['atom1'], bond['atom2']), 
                                    max(bond['atom1'], bond['atom2'])))
                print(f"⚠️ Removendo ligação excessiva do H{atom1['element']}")
            
            if atom2['element'] == 'H' and bond_count[bond['atom2']] > 1:
                bonds_to_remove.add((min(bond['atom1'], bond['atom2']), 
                                    max(bond['atom1'], bond['atom2'])))
                print(f"⚠️ Removendo ligação excessiva do H{atom2['element']}")
            
            # Regra 2: Verificar valência máxima
            val_max1 = valencias.get(atom1['element'], 4)
            val_max2 = valencias.get(atom2['element'], 4)
            
            if bond_count[bond['atom1']] > val_max1:
                bonds_to_remove.add((min(bond['atom1'], bond['atom2']), 
                                    max(bond['atom1'], bond['atom2'])))
                print(f"⚠️ Valência excedida para {atom1['element']}")
            
            if bond_count[bond['atom2']] > val_max2:
                bonds_to_remove.add((min(bond['atom1'], bond['atom2']), 
                                    max(bond['atom1'], bond['atom2'])))
                print(f"⚠️ Valência excedida para {atom2['element']}")
        
        # Filtrar ligações válidas
        bonds_validos = []
        for bond in bonds:
            key = (min(bond['atom1'], bond['atom2']), max(bond['atom1'], bond['atom2']))
            if key not in bonds_to_remove:
                bonds_validos.append(bond)
        
        if len(bonds_to_remove) > 0:
            print(f"✅ Removidas {len(bonds_to_remove)} ligações inválidas")
        
        return bonds_validos
    
    @staticmethod
    def auto_bonds(atoms):
        """Gera ligações baseado em distância atômica (valores otimizados)"""
        bonds = []
        
        # Raios covalentes aproximados (Angstroms)
        raios = {
            'H': 0.37, 'C': 0.77, 'N': 0.75, 'O': 0.73, 'F': 0.71,
            'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33,
            'Na': 1.54, 'Mg': 1.30, 'Ca': 1.74, 'Fe': 1.25, 'Cu': 1.28,
            'Zn': 1.22, 'Ag': 1.53, 'Au': 1.44, 'Hg': 1.49, 'Pb': 1.54
        }
        
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dx = atoms[i]['x'] - atoms[j]['x']
                dy = atoms[i]['y'] - atoms[j]['y']
                dz = atoms[i]['z'] - atoms[j]['z']
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                
                # Calcular distância máxima esperada baseada nos raios
                raio_i = raios.get(atoms[i]['element'], 0.8)
                raio_j = raios.get(atoms[j]['element'], 0.8)
                max_dist = (raio_i + raio_j) * 1.4  # 40% de margem
                min_dist = 0.3  # Distância mínima para evitar falsas ligações
                
                if min_dist < dist < max_dist:
                    # Determinar tipo baseado na distância
                    if dist < (raio_i + raio_j) * 1.1:
                        bond_type = 1  # Simples
                    elif dist < (raio_i + raio_j) * 1.2:
                        bond_type = 2  # Dupla
                    else:
                        bond_type = 1  # Simples
                    
                    bonds.append({'atom1': i, 'atom2': j, 'type': bond_type})
        
        # Validar ligações geradas automaticamente
        bonds = PubChemAPI.validar_ligacoes(atoms, bonds)
        
        return bonds
    
    @staticmethod
    def buscar_por_formula(formula):
        formula = formula.strip()
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{quote(formula)}/cids/txt"
        
        try:
            response = requests.get(url_cid, timeout=8)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                if lines and lines[0].isdigit():
                    cid = lines[0]
                    return PubChemAPI.buscar_por_cid(cid)
        except:
            pass
        return None
    
    @staticmethod
    def buscar_por_nome(nome):
        nome_traduzido = TradutorQuimico.traduzir(nome)
        
        variacoes = [
            nome_traduzido,
            nome_traduzido.lower(),
            nome_traduzido.capitalize(),
            nome_traduzido.replace(' ', ''),
        ]
        
        for variacao in set(variacoes):
            if len(variacao) < 2:
                continue
            
            try:
                url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(variacao)}/cids/txt"
                response = requests.get(url_cid, timeout=8)
                if response.status_code == 200:
                    lines = response.text.strip().split('\n')
                    if lines and lines[0].isdigit():
                        cid = lines[0]
                        return PubChemAPI.buscar_por_cid(cid)
            except:
                continue
        
        return None
    
    @staticmethod
    def buscar_por_cid(cid):
        url_sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        
        try:
            response = requests.get(url_sdf, headers={'Accept': 'chemical/x-mdl-sdfile'}, timeout=10)
            if response.status_code == 200:
                return PubChemAPI.parse_sdf(response.text)
        except:
            pass
        
        try:
            url_sdf_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
            response = requests.get(url_sdf_2d, headers={'Accept': 'chemical/x-mdl-sdfile'}, timeout=10)
            if response.status_code == 200:
                return PubChemAPI.parse_sdf(response.text)
        except:
            pass
        
        return None
    
    @staticmethod
    def buscar_automatico(query):
        query = query.strip()
        
        cached = cache_store.get(query)
        if cached:
            return cached
        
        resultado = None
        
        query_lower = query.lower()
        if query_lower in COMPOSTOS_ESPECIAIS:
            dados = COMPOSTOS_ESPECIAIS[query_lower]
            resultado = PubChemAPI.buscar_por_cid(dados['cid'])
            if resultado:
                cache_store.set(query, resultado)
                return resultado
        
        resultado = PubChemAPI.buscar_por_nome(query)
        if resultado:
            cache_store.set(query, resultado)
            return resultado
        
        resultado = PubChemAPI.buscar_por_formula(query)
        if resultado:
            cache_store.set(query, resultado)
            return resultado
        
        return None
    
    @staticmethod
    def parse_sdf(sdf_text):
        lines = sdf_text.strip().split('\n')
        if len(lines) < 4:
            return None
        
        atom_count, bond_count, start_line = 0, 0, 0
        
        # Procurar a linha de contagem
        for i in range(min(100, len(lines))):
            line = lines[i]
            if len(line) >= 6 and line[0:3].strip().isdigit():
                try:
                    atom_count = int(line[0:3].strip())
                    bond_count = int(line[3:6].strip())
                    start_line = i + 1
                    break
                except:
                    continue
        
        if atom_count == 0 or atom_count > 5000:
            return None
        
        # Pular linhas vazias
        while start_line < len(lines) and len(lines[start_line].strip()) < 5:
            start_line += 1
        
        # Parse dos átomos
        atoms = []
        atom_index_map = {}  # Mapear índice do SDF para índice da lista
        current_idx = 0
        
        for i in range(start_line, min(start_line + atom_count, len(lines))):
            line = lines[i]
            if len(line) >= 34:
                try:
                    x = float(line[0:10].strip())
                    y = float(line[10:20].strip())
                    z = float(line[20:30].strip())
                    element = line[31:34].strip()
                    
                    # Extrair elemento químico (ignorar números e caracteres especiais)
                    element_match = re.match(r'([A-Z][a-z]?)', element)
                    if element_match:
                        element = element_match.group(1)
                    
                    # Filtrar elementos válidos
                    if element and element[0].isalpha() and len(element) <= 2:
                        atom_index_map[current_idx] = len(atoms)
                        atoms.append({
                            'element': element, 
                            'x': float(x), 
                            'y': float(y), 
                            'z': float(z)
                        })
                    current_idx += 1
                except (ValueError, IndexError) as e:
                    print(f"Erro ao parsear átomo: {e}")
                    current_idx += 1
                    continue
        
        if len(atoms) == 0:
            print("Nenhum átomo válido encontrado no SDF")
            return None
        
        # Parse das ligações
        bonds = []
        bond_start = start_line + atom_count
        
        for i in range(bond_start, min(bond_start + bond_count, len(lines))):
            line = lines[i]
            if len(line) >= 9:
                try:
                    a1_original = int(line[0:3].strip()) - 1  # índices originais (1-based)
                    a2_original = int(line[3:6].strip()) - 1
                    bond_type = int(line[6:9].strip())
                    
                    # Mapear para índices atuais
                    if a1_original in atom_index_map and a2_original in atom_index_map:
                        a1 = atom_index_map[a1_original]
                        a2 = atom_index_map[a2_original]
                        
                        # Corrigir tipo de ligação (SDF às vezes tem 4=aromática, etc)
                        if bond_type == 4:
                            bond_type = 1  # Aromática trata como simples
                        elif bond_type > 3:
                            bond_type = 1
                        
                        if 0 <= a1 < len(atoms) and 0 <= a2 < len(atoms):
                            # Evitar duplicatas
                            if not any(b['atom1'] == a2 and b['atom2'] == a1 for b in bonds):
                                bonds.append({
                                    'atom1': a1, 
                                    'atom2': a2, 
                                    'type': bond_type
                                })
                except (ValueError, IndexError) as e:
                    print(f"Erro ao parsear ligação: {e}")
                    continue
        
        # IMPORTANTE: Validar ligações antes de retornar
        if bonds:
            print(f"🔗 {len(bonds)} ligações encontradas no SDF, validando...")
            bonds = PubChemAPI.validar_ligacoes(atoms, bonds)
            print(f"✅ {len(bonds)} ligações válidas após validação")
        elif len(atoms) > 1:
            print("⚠️ Nenhuma ligação encontrada, gerando automaticamente...")
            bonds = PubChemAPI.auto_bonds(atoms)
            print(f"✅ {len(bonds)} ligações geradas automaticamente")
        
        return {
            'atoms': atoms,
            'bonds': bonds,
            'atom_count': len(atoms),
            'bond_count': len(bonds)
        }

# ============================================
# ROTAS DA API
# ============================================

# ROTA PRINCIPAL - SERVE O HTML
@app.route('/')
def serve_index():
    return send_from_directory(os.path.dirname(__file__), 'index.html')

@app.route('/api/test', methods=['GET'])
def test():
    return jsonify({'status': 'ok', 'message': 'API Química na Vercel!'})

@app.route('/api/molecule/<path:query>', methods=['GET'])
def get_molecule(query):
    start_time = time.time()
    resultado = PubChemAPI.buscar_automatico(query)
    elapsed = time.time() - start_time
    
    if resultado:
        # Log para debug
        print(f"✅ Molécula encontrada: {query}")
        print(f"   - Átomos: {resultado['atom_count']}")
        print(f"   - Ligações: {resultado['bond_count']}")
        
        return jsonify({
            'success': True,
            'meta': {
                'query_original': query,
                'tempo_busca': f"{elapsed:.2f}s",
                'fonte': 'PubChem API + Dicionário PT'
            },
            'data': resultado
        })
    else:
        print(f"❌ Molécula não encontrada: {query}")
        return jsonify({
            'success': False,
            'error': f'Não foi possível encontrar: "{query}"'
        }), 404

@app.route('/api/traduzir/<path:texto>', methods=['GET'])
def test_traducao(texto):
    traduzido = TradutorQuimico.traduzir(texto)
    return jsonify({'original': texto, 'traduzido': traduzido})

@app.route('/api/compostos', methods=['GET'])
def listar_compostos_especiais():
    return jsonify({
        'total': len(COMPOSTOS_ESPECIAIS),
        'compostos': list(COMPOSTOS_ESPECIAIS.keys())
    })

@app.route('/api/cid/<int:cid>', methods=['GET'])
def get_by_cid(cid):
    resultado = PubChemAPI.buscar_por_cid(cid)
    if resultado:
        return jsonify({'success': True, 'data': resultado})
    return jsonify({'success': False, 'error': f'CID {cid} não encontrado'}), 404

# Handler para Vercel
app = app

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
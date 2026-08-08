/* ------------------------------------------------------------------
   meshio.js — read user-supplied 3D models into a plain triangle soup

   A wind tunnel needs geometry and nothing else: no materials, no UVs,
   no animation, no scene semantics. That lets us support the formats
   people actually download without pulling in a single dependency.

     .stl   binary + ASCII   — pure triangle soup
     .obj   text             — v / f, polygons fan-triangulated
     .glb   binary glTF 2.0  — geometry only, node transforms applied

   Every parser returns a Float32Array of length 9*ntri: three vertices
   per triangle, x,y,z each, in the file's own coordinate system.
------------------------------------------------------------------ */
(function (global) {
  'use strict';

  // ---- STL -----------------------------------------------------------
  function isBinarySTL(buf) {
    if (buf.byteLength < 84) return false;
    return 84 + new DataView(buf).getUint32(80, true) * 50 === buf.byteLength;
  }

  function parseSTLBinary(buf) {
    const dv = new DataView(buf);
    const n = dv.getUint32(80, true);
    const out = new Float32Array(n * 9);
    let o = 84, p = 0;
    for (let i = 0; i < n; i++) {
      o += 12;                                   // face normal — recomputed later
      for (let k = 0; k < 9; k++) { out[p++] = dv.getFloat32(o, true); o += 4; }
      o += 2;                                    // attribute byte count
    }
    return out;
  }

  function parseSTLAscii(txt) {
    const re = /vertex\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)/g;
    const nums = [];
    let m;
    while ((m = re.exec(txt))) nums.push(+m[1], +m[2], +m[3]);
    if (!nums.length) throw new Error('No vertices found in ASCII STL.');
    return new Float32Array(nums);
  }

  // ---- OBJ -----------------------------------------------------------
  function parseOBJ(txt) {
    const vx = [], vy = [], vz = [];
    const idx = [];
    const lines = txt.split('\n');
    for (let li = 0; li < lines.length; li++) {
      const line = lines[li].trim();
      if (line.length < 2) continue;
      const c0 = line[0];
      if (c0 === 'v' && line[1] === ' ') {
        const p = line.split(/\s+/);
        vx.push(+p[1]); vy.push(+p[2]); vz.push(+p[3]);
      } else if (c0 === 'f' && line[1] === ' ') {
        const p = line.split(/\s+/);
        const face = [];
        for (let i = 1; i < p.length; i++) {
          // "v", "v/vt", "v//vn", "v/vt/vn" — only the vertex index matters
          const v = parseInt(p[i], 10);
          if (!isNaN(v)) face.push(v < 0 ? vx.length + v : v - 1);
        }
        for (let i = 1; i + 1 < face.length; i++) {      // fan-triangulate
          idx.push(face[0], face[i], face[i + 1]);
        }
      }
    }
    if (!idx.length) throw new Error('No faces found in OBJ.');
    const out = new Float32Array(idx.length * 3);
    for (let i = 0; i < idx.length; i++) {
      const v = idx[i];
      out[i * 3] = vx[v]; out[i * 3 + 1] = vy[v]; out[i * 3 + 2] = vz[v];
    }
    return out;
  }

  // ---- glTF 2.0 binary ------------------------------------------------
  const CSIZE = { 5120: 1, 5121: 1, 5122: 2, 5123: 2, 5125: 4, 5126: 4 };
  const NCOMP = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4, MAT4: 16 };

  function readAccessor(g, bin, i) {
    const a = g.accessors[i];
    if (a.bufferView === undefined) return null;       // sparse-only: unsupported
    const bv = g.bufferViews[a.bufferView];
    const nc = NCOMP[a.type], cs = CSIZE[a.componentType];
    const stride = bv.byteStride || cs * nc;
    const base = bin.byteOffset + (bv.byteOffset || 0) + (a.byteOffset || 0);
    const dv = new DataView(bin.buffer);
    const out = a.componentType === 5126
      ? new Float32Array(a.count * nc) : new Uint32Array(a.count * nc);
    for (let k = 0; k < a.count; k++) {
      let p = base + k * stride;
      for (let c = 0; c < nc; c++) {
        let v;
        switch (a.componentType) {
          case 5126: v = dv.getFloat32(p, true); break;
          case 5125: v = dv.getUint32(p, true); break;
          case 5123: v = dv.getUint16(p, true); break;
          case 5122: v = dv.getInt16(p, true); break;
          case 5121: v = dv.getUint8(p); break;
          default: v = dv.getInt8(p);
        }
        out[k * nc + c] = v;
        p += cs;
      }
    }
    return out;
  }

  const IDENT = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];

  function mul4(a, b) {                    // column-major, as glTF stores them
    const o = new Array(16);
    for (let c = 0; c < 4; c++) {
      for (let r = 0; r < 4; r++) {
        o[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] +
          a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3];
      }
    }
    return o;
  }

  function nodeMatrix(n) {
    if (n.matrix) return n.matrix.slice();
    const t = n.translation || [0, 0, 0];
    const r = n.rotation || [0, 0, 0, 1];          // quaternion x,y,z,w
    const s = n.scale || [1, 1, 1];
    const [x, y, z, w] = r;
    const x2 = x + x, y2 = y + y, z2 = z + z;
    const xx = x * x2, xy = x * y2, xz = x * z2;
    const yy = y * y2, yz = y * z2, zz = z * z2;
    const wx = w * x2, wy = w * y2, wz = w * z2;
    return [
      (1 - (yy + zz)) * s[0], (xy + wz) * s[0], (xz - wy) * s[0], 0,
      (xy - wz) * s[1], (1 - (xx + zz)) * s[1], (yz + wx) * s[1], 0,
      (xz + wy) * s[2], (yz - wx) * s[2], (1 - (xx + yy)) * s[2], 0,
      t[0], t[1], t[2], 1
    ];
  }

  function parseGLB(buf) {
    const dv = new DataView(buf);
    if (dv.getUint32(0, true) !== 0x46546C67) throw new Error('Not a GLB file (bad magic).');
    const total = Math.min(dv.getUint32(8, true), buf.byteLength);
    let o = 12, g = null, bin = null;
    while (o + 8 <= total) {
      const len = dv.getUint32(o, true), type = dv.getUint32(o + 4, true);
      const start = o + 8;
      if (type === 0x4E4F534A) {
        g = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, start, len)));
      } else if (type === 0x004E4942) {
        bin = new Uint8Array(buf, start, len);
      }
      o = start + len + ((4 - (len & 3)) & 3);
    }
    if (!g) throw new Error('GLB contains no JSON chunk.');
    if (!bin) throw new Error('GLB has no binary chunk (external .bin files are not supported).');
    if (g.extensionsRequired && g.extensionsRequired.length) {
      throw new Error('This GLB needs ' + g.extensionsRequired.join(', ') +
        '. Re-export without compression (uncheck Draco / meshopt) and try again.');
    }
    if (!g.meshes || !g.meshes.length) throw new Error('GLB contains no meshes.');

    const chunks = [];
    let count = 0;

    const emit = (prim, M) => {
      if (prim.mode !== undefined && prim.mode !== 4) return;      // triangles only
      const pa = prim.attributes && prim.attributes.POSITION;
      if (pa === undefined) return;
      const pos = readAccessor(g, bin, pa);
      if (!pos) return;
      let ind = prim.indices !== undefined ? readAccessor(g, bin, prim.indices) : null;
      const ntri = ind ? ind.length / 3 : pos.length / 9;
      const tri = new Float32Array(ntri * 9);
      for (let t = 0; t < ntri * 3; t++) {
        const v = (ind ? ind[t] : t) * 3;
        const x = pos[v], y = pos[v + 1], z = pos[v + 2];
        tri[t * 3] = M[0] * x + M[4] * y + M[8] * z + M[12];
        tri[t * 3 + 1] = M[1] * x + M[5] * y + M[9] * z + M[13];
        tri[t * 3 + 2] = M[2] * x + M[6] * y + M[10] * z + M[14];
      }
      chunks.push(tri); count += tri.length;
    };

    const walk = (ni, parent, depth) => {
      if (depth > 64) return;
      const n = g.nodes[ni];
      if (!n) return;
      const M = mul4(parent, nodeMatrix(n));
      if (n.mesh !== undefined && g.meshes[n.mesh]) {
        for (const prim of g.meshes[n.mesh].primitives) emit(prim, M);
      }
      if (n.children) for (const c of n.children) walk(c, M, depth + 1);
    };

    const scene = g.scenes && g.scenes[g.scene === undefined ? 0 : g.scene];
    if (scene && scene.nodes) {
      for (const ni of scene.nodes) walk(ni, IDENT, 0);
    } else if (g.nodes) {
      for (let i = 0; i < g.nodes.length; i++) walk(i, IDENT, 0);
    } else {
      for (const m of g.meshes) for (const prim of m.primitives) emit(prim, IDENT);
    }

    if (!count) throw new Error('GLB has meshes but no readable triangles.');
    const out = new Float32Array(count);
    let p = 0;
    for (const c of chunks) { out.set(c, p); p += c.length; }
    return out;
  }

  // ---- attribution -----------------------------------------------------
  // Most model sites ship credits inside the file. glTF has a standard place
  // for it (asset.copyright and asset.extras, which is what Sketchfab fills
  // in); OBJ and STL only have comment headers, so those are best-effort.
  //
  // Everything here is untrusted text from a downloaded file, so the caller
  // must render it as text and validate any URL before making it a link.

  function pickFields(txt) {
    const m = {};
    const grab = (key) => {
      const re = new RegExp('^[#\\s]*' + key + '\\s*[:=]\\s*(.+)$', 'im');
      const hit = re.exec(txt);
      return hit ? hit[1].trim() : undefined;
    };
    m.title = grab('title') || grab('name');
    m.author = grab('author') || grab('creator') || grab('artist') || grab('by');
    m.license = grab('licen[cs]e');
    m.source = grab('source') || grab('url');
    return m;
  }

  function metaFromGLTF(g) {
    const a = (g && g.asset) || {};
    const x = a.extras || {};
    const m = {
      title: x.title, author: x.author, license: x.license,
      source: x.source, copyright: a.copyright, generator: a.generator
    };
    // some exporters put everything in one copyright string instead
    if (a.copyright && !m.author && !m.license) {
      Object.assign(m, pickFields(a.copyright));
    }
    return m;
  }

  /** Read just the JSON chunk of a GLB — no geometry decoding. */
  function glbJSON(buf) {
    const dv = new DataView(buf);
    if (buf.byteLength < 20 || dv.getUint32(0, true) !== 0x46546C67) return null;
    const total = Math.min(dv.getUint32(8, true), buf.byteLength);
    let o = 12;
    while (o + 8 <= total) {
      const len = dv.getUint32(o, true), type = dv.getUint32(o + 4, true);
      if (type === 0x4E4F534A) {
        try {
          return JSON.parse(new TextDecoder().decode(new Uint8Array(buf, o + 8, len)));
        } catch (e) { return null; }
      }
      o = o + 8 + len + ((4 - (len & 3)) & 3);
    }
    return null;
  }

  const printable = (s) => s.replace(/[^\x20-\x7E -￿]+/g, ' ').trim();

  /**
   * Best-effort credits for a model file.
   * @returns {title,author,license,source,copyright,generator,notes} — any
   *          field may be missing; returns null when the file says nothing.
   */
  function meta(name, buf) {
    const ext = (name.split('.').pop() || '').toLowerCase();
    let m = null;

    if (ext === 'glb' || (buf.byteLength > 12 &&
      new DataView(buf).getUint32(0, true) === 0x46546C67)) {
      const g = glbJSON(buf);
      if (g) m = metaFromGLTF(g);
    } else if (ext === 'gltf') {
      try { m = metaFromGLTF(JSON.parse(new TextDecoder().decode(new Uint8Array(buf)))); }
      catch (e) { /* not our problem here */ }
    } else if (ext === 'stl' && isBinarySTL(buf)) {
      const head = printable(new TextDecoder().decode(new Uint8Array(buf, 0, 80)));
      if (head) m = Object.assign({ notes: head }, pickFields(head));
    } else {
      // OBJ (and ASCII STL): leading comment block
      const head = new TextDecoder().decode(
        new Uint8Array(buf, 0, Math.min(4096, buf.byteLength)));
      const lines = [];
      for (const raw of head.split('\n')) {
        const l = raw.trim();
        if (l.startsWith('#')) { lines.push(printable(l.replace(/^#+\s?/, ''))); }
        else if (l) break;                    // comments only run at the top
        if (lines.length > 12) break;
      }
      const block = lines.filter(Boolean).join('\n');
      if (block) m = Object.assign({ notes: block }, pickFields(block));
    }

    if (!m) return null;
    for (const k of Object.keys(m)) {
      if (typeof m[k] === 'string') {
        m[k] = printable(m[k]).slice(0, 300);
        if (!m[k]) delete m[k];
      } else if (m[k] === undefined) delete m[k];
    }
    return Object.keys(m).length ? m : null;
  }

  // ---- dispatch --------------------------------------------------------
  /**
   * @param name  file name (used for the extension)
   * @param buf   ArrayBuffer of the whole file
   * @returns Float32Array, 9 floats per triangle
   */
  function parse(name, buf) {
    const ext = (name.split('.').pop() || '').toLowerCase();
    if (ext === 'glb') return parseGLB(buf);
    if (ext === 'gltf') {
      throw new Error('Plain .gltf references external files. Export as .glb (binary glTF) instead.');
    }
    if (ext === 'stl') {
      return isBinarySTL(buf) ? parseSTLBinary(buf)
        : parseSTLAscii(new TextDecoder().decode(new Uint8Array(buf)));
    }
    if (ext === 'obj') return parseOBJ(new TextDecoder().decode(new Uint8Array(buf)));
    // no extension match — sniff
    if (buf.byteLength > 12 && new DataView(buf).getUint32(0, true) === 0x46546C67) return parseGLB(buf);
    if (isBinarySTL(buf)) return parseSTLBinary(buf);
    const txt = new TextDecoder().decode(new Uint8Array(buf, 0, Math.min(4096, buf.byteLength)));
    if (/^\s*solid/.test(txt)) return parseSTLAscii(new TextDecoder().decode(new Uint8Array(buf)));
    if (/^\s*(v|vn|vt|f|o|g|mtllib)\s/m.test(txt)) return parseOBJ(new TextDecoder().decode(new Uint8Array(buf)));
    throw new Error('Unrecognised file. Supported: .obj, .stl, .glb');
  }

  global.MeshIO = {
    parse, meta, parseOBJ, parseGLB, parseSTLBinary, parseSTLAscii, isBinarySTL
  };
})(window);

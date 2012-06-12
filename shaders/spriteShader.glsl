// lightingTexturingShader.glsl
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com



#ifdef _VERTEX_

// Attributes : per vertex data
in vec4 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexUvs;
in vec4 vertexColor;

// Varyings : data to transmit to fragments
smooth out vec4 position;
smooth out vec4 normal;
smooth out vec2 uvs;
smooth out vec4 localColor;

void main()
{
    if (filledData[0]) position = model * vertexPosition;
    if (filledData[1]) normal = model * vec4(vertexNormal, 0.0);
    if (filledData[2]) uvs = vertexUvs;
    if (filledData[3]) localColor = vertexColor;

    gl_PointSize=15;

    gl_Position = projection * view * model * vertexPosition;
}


#endif




#ifdef _FRAGMENT_

// Varyings : data receved and interpolated from the vertex shaders
smooth in vec4 position;
smooth in vec4 normal;
smooth in vec2 uvs;
smooth in vec4 localColor;

// Final output
out vec4 fragColor;

void main()
{
    vec3 text=texture(textureUnitDiffuse, gl_PointCoord).rgb;
    float alpha=localColor.a*((text.r+text.g+text.b)/3.0);
    fragColor=vec4(localColor.rgb*alpha, alpha);
}


#endif

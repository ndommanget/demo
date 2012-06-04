// lightingShader.glsl
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
    vec4 diffuseColorMix=vec4(material.diffuse);
    // If color
    if (filledData[3]) diffuseColorMix=localColor;
    // diffuseColorMix=mix(localColor, material.diffuse, 0.3);
    
    // If no normal
    if (!filledData[1]) fragColor = diffuseColorMix;
    else
    {
        vec4 L=normalize(light.position); // Direction of light from fragment -> light.position[3]==0.0 : Directional light
        if (light.position.w==1.0) L=normalize(light.position-position); //   -> light.position[3]==1.0 : Point light
        vec4 V=normalize(center-position); // Direction from fragment to camera center
        vec4 R=normalize(reflect(-L, normal)); // Direction of reflected light beam, from fragment
        vec4 N=normalize(normal); // Normal
    
        float ambientValue=light.power;
        float diffuseValue=light.power * max( dot(L, N), 0.0);
        float specularValue=light.power * pow( max(dot(R, V), 0.0), material.shininess);
        vec4 ambientContribution=vec4(ambientValue*material.ka*material.ambient.rgb, material.ambient.a);
        vec4 diffuseContribution=vec4(diffuseValue*material.kd*diffuseColorMix.rgb, diffuseColorMix.a);
        vec4 specularContribution=vec4(specularValue*material.ks*material.specular.rgb, material.specular.a);
  
        fragColor = ambientContribution + diffuseContribution + specularContribution;
    }
}


#endif
